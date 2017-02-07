/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: CMap2.hpp
 *
 * Description:
 *      description
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2017-01-31   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_CMAP2_HPP_
#define SRC_CMAP2_HPP_

#include "CMAPLib.hpp"
#include "query_parser.hpp"

#include "cpplinq.hpp"
#include "likely.h"

#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstddef>

#undef NSIG

using size_type = std::size_t;

#include <sys/time.h>
long int timestamp()
{
    timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec;// + tv.tv_usec / 1e6;
}

struct CMAP2Updated
{
    enum {NGENES = 10174};

    CMAP2Updated(size_type nsig = 476251, size_type nskip = 0, std::string path = "", std::vector<double> const * gt = nullptr) :
        NSIG(nsig), NSKIP(nskip), PATH(path), m_gt(gt)
    {}

    int init(std::vector<int> const & genes);
    std::vector<double> getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);
    std::vector<double> getWTKScomb_(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);

    size_type const NSIG;
    size_type const NSKIP;
    std::string const PATH;
    std::vector<double> const * m_gt;
    std::unordered_map<int, int> m_gene_to_idx;
    CMAPLib m_cmap_lib;
};

int
CMAP2Updated::init(std::vector<int> const & genes)
{
    m_gene_to_idx.clear();
    namespace q = ::cpplinq;
    q::from(genes)
        >> q::zip_with(q::detail::int_range(0, INT_MAX))
        >> q::aggregate(&m_gene_to_idx,
            [](std::unordered_map<int, int> * seed, std::pair<int, int> const & kv)
            {
                seed->emplace(kv.first, kv.second);
                return seed;
            })
        ;

    // cache warm up ??? TODO

    return 0;
}


std::vector<double>
CMAP2Updated::getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn)
{
    std::cerr << "Will process " << NSIG - NSKIP << " signatures, [" << NSKIP << ',' << NSIG << ")\n";

    const auto time0 = timestamp();

    assert(q_up.size() == q_dn.size());
    auto const NQRY = q_up.size();

    auto q_indexer = [this](std::vector<std::string> & qry)
        {
            namespace q = ::cpplinq;
            return q::from(qry)
                >> q::select([this](std::string const & s)
                    {
                        auto q_parsed = parse_query(s);
                        std::transform(q_parsed.begin(), q_parsed.end(), q_parsed.begin(),
                            [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });
                        return q_parsed;
                    })
                >> q::to_vector();
        };

    auto const q_up_indexed = q_indexer(q_up);
    auto const q_dn_indexed = q_indexer(q_dn);

    std::vector<double> sig;
    std::vector<int> ranks;

    std::vector<double> ret;
    ret.reserve(q_up.size() * NSIG);

    for (auto six = NSKIP; six < NSIG; ++six)
    {
        sig = m_cmap_lib.loadFromDoubleFile(PATH + "scoresBySig", six * NGENES, NGENES);
        ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySig", six * NGENES, NGENES);
        // cache i/o ???

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto & q_up_indices = q_up_indexed[qix];
            auto & q_dn_indices = q_dn_indexed[qix];

            auto const QUSIZE = q_up_indices.size();
            auto const QDSIZE = q_dn_indices.size();

            std::vector<std::pair<std::uint16_t, std::uint16_t>> up_ranks;
            up_ranks.reserve(QUSIZE);
            double Sum_up_abs_scores = 0.;
            std::transform(q_up_indices.cbegin(), q_up_indices.cend(), std::back_inserter(up_ranks),
                [&ranks, &sig, &Sum_up_abs_scores](int ix)
                {
                    auto abs_score = std::abs(sig[ix]);
                    Sum_up_abs_scores += abs_score;
                    return std::make_pair<std::uint16_t, std::uint16_t>(ranks[ix], ix);
                });

            std::sort(up_ranks.begin(), up_ranks.end(),
                [](std::pair<std::uint16_t, std::uint16_t> const & p, std::pair<std::uint16_t, std::uint16_t> const & q)
                { return p.first < q.first;});

            double const upenalty = -1. / (NGENES - QUSIZE);
            double wkts_up = 0.;
            {
                double _acc = 0.;
                double _extr = 0.;
                int prev = 0; // rank indices are 1-based
                for (auto const & rs : up_ranks)
                {
                    _acc += (rs.first - prev - 1) * upenalty;
                    auto abs_acc = std::abs(_acc);
                    if (UNLIKELY(abs_acc > _extr))
                    {
                        _extr = abs_acc;
                        wkts_up = _acc;
                    }
                    prev = rs.first;

                    _acc += std::abs(sig[rs.second]) / Sum_up_abs_scores;
                    abs_acc = std::abs(_acc);
                    if (UNLIKELY(abs_acc > _extr))
                    {
                        _extr = abs_acc;
                        wkts_up = _acc;
                    }
                }
            }

            std::vector<std::pair<std::uint16_t, std::uint16_t>> dn_ranks;
            dn_ranks.reserve(QDSIZE);
            double Sum_dn_abs_scores = 0.;
            std::transform(q_dn_indices.cbegin(), q_dn_indices.cend(), std::back_inserter(dn_ranks),
                [&ranks, &sig, &Sum_dn_abs_scores](int ix)
                {
                    auto abs_score = std::abs(sig[ix]);
                    Sum_dn_abs_scores += abs_score;
                    return std::make_pair<std::uint16_t, std::uint16_t>(ranks[ix], ix);
                });

            std::sort(dn_ranks.begin(), dn_ranks.end(),
                [](std::pair<std::uint16_t, std::uint16_t> const & p, std::pair<std::uint16_t, std::uint16_t> const & q)
                { return p.first < q.first;});

            double const dpenalty = -1. / (NGENES - QDSIZE);
            double wkts_dn = 0.;
            {
                double _acc = 0.;
                double _extr = 0.;
                int prev = 0; // rank indices are 1-based
                for (auto const & rs : dn_ranks)
                {
                    _acc += (rs.first - prev - 1) * dpenalty;
                    auto abs_acc = std::abs(_acc);
                    if (UNLIKELY(abs_acc > _extr))
                    {
                        _extr = abs_acc;
                        wkts_dn = _acc;
                    }
                    prev = rs.first;

                    _acc += std::abs(sig[rs.second]) / Sum_dn_abs_scores;
                    abs_acc = std::abs(_acc);
                    if (UNLIKELY(abs_acc > _extr))
                    {
                        _extr = abs_acc;
                        wkts_dn = _acc;
                    }
                }
            }

            auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;
            //std::cout << "wkts : " << wkts << "\n";
            ret.push_back(wkts);
            if (UNLIKELY(m_gt != nullptr))
            {
                auto ref = (*m_gt)[(six - NSKIP) * NQRY + qix];
                if (std::abs(ref - wkts) >= 0.001)
                {
                    std::cout << "Sig: " << six + 1 << ", Query: " << qix + 1 \
                        << " Ref: " << ref << " WKTS: " << wkts << " up/dn " << wkts_up << ' ' << wkts_dn << '\n';
                }
            }
        }
    }

    return ret;
}

#endif /* SRC_CMAP2_HPP_ */
