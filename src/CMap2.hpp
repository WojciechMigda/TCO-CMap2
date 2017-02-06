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

    CMAP2Updated(size_type nsig = 476251, std::string path = "", std::vector<double> const * gt = nullptr) :
        NSIG(nsig), PATH(path), m_gt(gt)
    {}

    int init(std::vector<int> const & genes);
    std::vector<double> getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);
    std::vector<double> getWTKScomb_(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);

    size_type const NSIG;
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
    const auto time0 = timestamp();

    assert(q_up.size() == q_dn.size());
    auto const NQRY = q_up.size();

    std::vector<bool> qmasks[2 * NQRY];
    for (auto & m : qmasks)
    {
        m.resize(NGENES);
    }

    std::vector<std::vector<std::uint32_t>> q_up_indices(NQRY);
    std::vector<std::vector<std::uint32_t>> q_dn_indices(NQRY);

    for (auto qix = 0u; qix < NQRY; ++qix)
    {
        q_up_indices[qix] = parse_query(q_up[qix]);
        auto & q_up_parsed = q_up_indices[qix];
        std::transform(q_up_parsed.begin(), q_up_parsed.end(), q_up_parsed.begin(),
            [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });

        q_dn_indices[qix] = parse_query(q_dn[qix]);
        auto & q_dn_parsed = q_dn_indices[qix];
        std::transform(q_dn_parsed.begin(), q_dn_parsed.end(), q_dn_parsed.begin(),
            [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });

        for (auto gix : q_up_parsed)
        {
            qmasks[2 * qix + 0][gix] = true;
        }
        for (auto gix : q_dn_parsed)
        {
            qmasks[2 * qix + 1][gix] = true;
        }
    }

    std::vector<double> ret;
    ret.reserve(q_up.size() * NSIG);

    std::vector<double> sig;
    std::vector<int> ranks;

    for (auto six = 0u; six < NSIG; ++six)
    {
        sig = m_cmap_lib.loadFromDoubleFile(PATH + "scoresBySig", six * NGENES, NGENES);
        //ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySig", six * NGENES, NGENES);

        namespace q = ::cpplinq;
        auto scores_ranks = q::from(sig)
            //>> q::select([](double x){ return std::abs(x); })
            >> q::zip_with(q::detail::int_range(0, INT_MAX))
            >> q::orderby_descending([](std::pair<double, int> const & p){ return /*std::abs*/(p.first); })
            //>> q::select([](std::pair<double, int> const & p){ return p.second; })
            >> q::to_vector(sig.size())
            ;


        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto const Sum_up_abs_scores = std::accumulate(
                q_up_indices[qix].cbegin(), q_up_indices[qix].cend(), 0.,
                [&sig](double acc, int ix){ return acc + std::abs(sig[ix]); });
            auto const Sum_dn_abs_scores = std::accumulate(
                q_dn_indices[qix].cbegin(), q_dn_indices[qix].cend(), 0.,
                [&sig](double acc, int ix){ return acc + std::abs(sig[ix]); });

            auto const QUSIZE = q_up_indices[qix].size();
            auto const QDSIZE = q_dn_indices[qix].size();

            double const upenalty = -1. / (NGENES - QUSIZE);
            double const dpenalty = -1. / (NGENES - QDSIZE);

            double wkts_up = 0.;
            double acc_up = 0.;
            double extr_up = 0.;

            double wkts_dn = 0.;
            double acc_dn = 0.;
            double extr_dn = 0.;

            for (auto ix = 0u; ix < NGENES; ++ix)
            {
                auto const & score_rank = scores_ranks[ix];
                auto const gix = score_rank.second;

                bool is_up_qry = qmasks[qix * 2 + 0][gix];
                bool is_dn_qry = qmasks[qix * 2 + 1][gix];

                if (is_up_qry)
                {
                    acc_up += std::abs(score_rank.first) / Sum_up_abs_scores;
                }
                else
                {
                    acc_up += upenalty;
                }

                if (is_dn_qry)
                {
                    acc_dn += std::abs(score_rank.first) / Sum_dn_abs_scores;
                }
                else
                {
                    acc_dn += dpenalty;
                }

                if (std::abs(acc_up) > extr_up)
                {
                    extr_up = std::abs(acc_up);
                    wkts_up = acc_up;
                }
                if (std::abs(acc_dn) > extr_dn)
                {
                    extr_dn = std::abs(acc_dn);
                    wkts_dn = acc_dn;
                }
            }
            auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;

            ret.push_back(wkts);
            if (m_gt)
            {
                auto ref = (*m_gt)[six * NQRY + qix];
                if (std::abs(ref - wkts) >= 0.001)
                {
                    std::cout << "Ref: " << ref << " WKTS: " << wkts << " up/dn " << wkts_up << ' ' << wkts_dn << '\n';
                }
            }
        }
    }
    return ret;
}

#endif /* SRC_CMAP2_HPP_ */
