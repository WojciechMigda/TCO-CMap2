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
#include <cmath>

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

    double runs[NGENES];

    for (auto six = 0u; six < NSIG; ++six)
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


            {
                double const upenalty = -1. / (NGENES - QUSIZE);
                double Sum_up_abs_scores = 0.;

                std::fill(std::begin(runs), std::end(runs), NAN);
                for (auto const ix : q_up_indices)
                {
                    auto abs_score = std::abs(sig[ix]);
                    runs[ranks[ix] - 1] = abs_score;
                    Sum_up_abs_scores += abs_score;
                }

                double wkts_up = 0.;
                double _acc = 0.;
                double _extr = 0.;
                for (auto const s : runs)
                {
                    if (std::isnan(s))
                    {
                        _acc += upenalty;
                    }
                    else
                    {
                        _acc += s / Sum_up_abs_scores;
                    }
                    auto abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_up = _acc;
                    }
                }

                double const dpenalty = -1. / (NGENES - QDSIZE);
                double Sum_dn_abs_scores = 0.;

                std::fill(std::begin(runs), std::end(runs), NAN);
                for (auto const ix : q_dn_indices)
                {
                    auto abs_score = std::abs(sig[ix]);
                    runs[ranks[ix] - 1] = abs_score;
                    Sum_dn_abs_scores += abs_score;
                }

                double wkts_dn = 0.;
                _acc = 0.;
                _extr = 0.;
                for (auto const s : runs)
                {
                    if (std::isnan(s))
                    {
                        _acc += dpenalty;
                    }
                    else
                    {
                        _acc += s / Sum_dn_abs_scores;
                    }
                    auto abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_dn = _acc;
                    }
                }

                auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;
                //std::cout << "wkts : " << wkts << "\n";
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

//            std::vector<std::pair<int, double>> up_ranks;
//            up_ranks.reserve(QUSIZE);
//            double Sum_up_abs_scores = 0.;
//            std::transform(q_up_indices.cbegin(), q_up_indices.cend(), std::back_inserter(up_ranks),
//                [&ranks, &sig, &Sum_up_abs_scores](int ix)
//                {
//                    auto abs_score = std::abs(sig.at(ix));
//                    Sum_up_abs_scores += abs_score;
//                    return std::make_pair(ranks.at(ix), abs_score);
//                });
//
//            std::vector<std::pair<int, double>> dn_ranks;
//            dn_ranks.reserve(QDSIZE);
//            double Sum_dn_abs_scores = 0.;
//            std::transform(q_dn_indices.cbegin(), q_dn_indices.cend(), std::back_inserter(dn_ranks),
//                [&ranks, &sig, &Sum_dn_abs_scores](int ix)
//                {
//                    auto abs_score = std::abs(sig.at(ix));
//                    Sum_dn_abs_scores += abs_score;
//                    return std::make_pair(ranks.at(ix), abs_score);
//                });
//
//            std::qsort(&up_ranks[0], up_ranks.size(), sizeof (up_ranks.front()),
//                [](const void * avp, const void * bvp)
//                {
//                    auto ap = static_cast<std::pair<int, double> const *>(avp);
//                    auto bp = static_cast<std::pair<int, double> const *>(bvp);
//                    return (ap->first > bp->first) - (ap->first < bp->first);
//                }
//            );
//            std::qsort(&dn_ranks[0], dn_ranks.size(), sizeof (dn_ranks.front()),
//                [](const void * avp, const void * bvp)
//                {
//                    auto ap = static_cast<std::pair<int, double> const *>(avp);
//                    auto bp = static_cast<std::pair<int, double> const *>(bvp);
//                    return (ap->first > bp->first) - (ap->first < bp->first);
//                }
//            );
//
//            double const upenalty = -1. / (NGENES - QUSIZE);
//            double wkts_up = 0.;
//            {
//                double _acc = 0.;
//                double _extr = 0.;
//                int prev = 0; // rank indices are 1-based
//                for (auto const & rs : up_ranks)
//                {
//                    _acc += (rs.first - prev - 1) * upenalty;
//                    auto abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_up = _acc;
//                    }
//                    prev = rs.first;
//
//                    _acc += rs.second / Sum_up_abs_scores;
//                    abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_up = _acc;
//                    }
//                }
//            }
//
//            double const dpenalty = -1. / (NGENES - QDSIZE);
//            double wkts_dn = 0.;
//            {
//                double _acc = 0.;
//                double _extr = 0.;
//                int prev = 0; // rank indices are 1-based
//                for (auto const & rs : dn_ranks)
//                {
//                    _acc += (rs.first - prev - 1) * dpenalty;
//                    auto abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_dn = _acc;
//                    }
//                    prev = rs.first;
//
//                    _acc += rs.second / Sum_dn_abs_scores;
//                    abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_dn = _acc;
//                    }
//                }
//            }
//
//            auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;
//            //std::cout << "wkts : " << wkts << "\n";
//            ret.push_back(wkts);
//            if (m_gt)
//            {
//                auto ref = (*m_gt)[six * NQRY + qix];
//                if (std::abs(ref - wkts) >= 0.001)
//                {
//                    std::cout << "Ref: " << ref << " WKTS: " << wkts << " up/dn " << wkts_up << ' ' << wkts_dn << '\n';
//                }
//            }
        }
    }

    return ret;
}

#endif /* SRC_CMAP2_HPP_ */
