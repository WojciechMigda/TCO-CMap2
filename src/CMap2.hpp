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
#include <functional>

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
    enum {RANK_PER_BUCKET = 1024};
    enum {NBUCKET = (NGENES + RANK_PER_BUCKET - 1) / RANK_PER_BUCKET};

    CMAP2Updated(size_type nsig = 476251, size_type nskip = 0, std::string path = "", std::vector<double> const * gt = nullptr) :
        NSIG(nsig), NSKIP(nskip), PATH(path), m_gt(gt)
    {}

    int init(std::vector<int> const & genes);
    std::vector<double> getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);
    std::vector<double> getWTKScomb_old(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);

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

    assert(q_up.size() == q_dn.size());
    auto const NQRY = q_up.size();

    std::vector<double> ret;
    ret.reserve(NQRY * NSIG);

    using gene_index_t = std::uint16_t;
    using score_index_t = std::uint16_t;

    using query_indexed_t = std::vector<gene_index_t>;

    auto q_indexer = [this](std::vector<std::string> & qry) -> std::vector<query_indexed_t>
        {
            namespace q = ::cpplinq;
            return q::from(qry)
                >> q::select([this](std::string const & s)
                    {
                        auto q_parsed = parse_query(s);
                        std::transform(q_parsed.begin(), q_parsed.end(), q_parsed.begin(),
                            [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });
                        return query_indexed_t(q_parsed.cbegin(), q_parsed.cend());
                    })
                >> q::to_vector();
        };

    // indices to unsorted scores
    auto const q_up_indexed = q_indexer(q_up);
    auto const q_dn_indexed = q_indexer(q_dn);

    typedef struct
    {
        std::vector<std::pair<double, score_index_t>> scores;
        double sum;
    } query_state_t;

    std::vector<query_state_t> q_up_states(NQRY);
    std::vector<query_state_t> q_dn_states(NQRY);

    std::vector<query_state_t *> gene_buckets[NGENES];

    for (auto qix = 0u; qix < NQRY; ++qix)
    {
        for (auto const gene_ix : q_up_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_up_states[qix]);
        }
        for (auto const gene_ix : q_dn_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_dn_states[qix]);
        }
    }

    for (auto six = 0u; six < NSIG; ++six)
    {
        auto sigs = m_cmap_lib.loadFromDoubleFile(PATH + "scoresBySig", six * NGENES, NGENES);

//        namespace q = ::cpplinq;
//        auto ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySig", six * NGENES, NGENES);
//        auto inv_ranks =
//            q::from(ranks)
//            >> q::zip_with(q::detail::int_range(0, INT_MAX))
//            >> q::orderby_ascending([](std::pair<int, int> const & ri){ return ri.first; })
//            >> q::select([](std::pair<int, int> const & ri){ return ri.second; })
//            >> q::to_vector(NGENES)
//            ;
        auto _ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySigInv", six * NGENES / 2, NGENES / 2);
        auto inv_ranks = (std::uint16_t const *)_ranks.data();

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            q_up_states[qix].scores.clear();
            q_up_states[qix].sum = 0.;

            q_dn_states[qix].scores.clear();
            q_dn_states[qix].sum = 0.;
        }

        for (auto gix = 0; gix < NGENES; ++gix)
        {
            auto inv_rank = inv_ranks[gix];
            auto score = std::abs(sigs[inv_rank]);

            for (auto q_state_p : gene_buckets[inv_rank])
            {
                q_state_p->scores.emplace_back(score, gix);
                q_state_p->sum += score;
            }
        }

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto wtks_calc = [](query_state_t const & q_state)
                {
                    auto QSIZE = q_state.scores.size();
                    double penalty = -1. / (NGENES - QSIZE);

                    double wtks = 0.;
                    double _acc = 0.;
                    double _min = 0.;
                    double _max = 0.;
                    int prev = 0; // rank indices are 1-based

                    for (auto const & sr : q_state.scores)
                    {
                        _acc += (sr.second - prev - 0) * penalty;
                        _min = std::min(_min, _acc);

                        prev = sr.second + 1;

                        _acc += sr.first / q_state.sum;
                        _max = std::max(_max, _acc);
                    }
                    if (_max > std::abs(_min))
                    {
                        wtks = _max;
                    }
                    else
                    {
                        wtks = _min;
                    }
                    return wtks;
                };

            double wtks_up = wtks_calc(q_up_states[qix]);
            double wtks_dn = wtks_calc(q_dn_states[qix]);

            auto wtks = (wtks_dn * wtks_up < 0.) ? (wtks_up - wtks_dn) / 2 : 0.;
            ret.push_back(wtks);

            if (UNLIKELY(m_gt != nullptr))
            {
                auto ref = (*m_gt)[(six - NSKIP) * NQRY + qix];
                if (std::abs(ref - wtks) >= 0.001)
                {
                    std::cout << "Sig: " << six + 1 << ", Query: " << qix + 1 \
                        << " Ref: " << ref << " WTKS: " << wtks << " up/dn " << wtks_up << ' ' << wtks_dn << '\n';
                }
            }
        }
    }

    return ret;
}

std::vector<double>
CMAP2Updated::getWTKScomb_old(std::vector<std::string> & q_up, std::vector<std::string> & q_dn)
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
                        return std::vector<std::uint16_t>(q_parsed.cbegin(), q_parsed.cend());
                        //return q_parsed;
                    })
                >> q::to_vector();
        };

    auto const q_up_indexed = q_indexer(q_up);
    auto const q_dn_indexed = q_indexer(q_dn);

    std::vector<double> sig;
    std::vector<int> ranks;

    std::vector<double> ret;
    ret.reserve(NQRY * NSIG);

    typedef union
    {
        struct
        {
            std::uint16_t score_ix;
            std::uint16_t rank;
        } rs;
        std::uint32_t packed;
    } rank_score_t;

    std::vector<rank_score_t> hranks[NBUCKET];
    for (auto & bucket : hranks)
    {
        bucket.reserve(250 / (2 * NBUCKET));
    }

    for (auto six = NSKIP; six < NSIG; ++six)
    {
        sig = m_cmap_lib.loadFromDoubleFile(PATH + "scoresBySig", six * NGENES, NGENES);
        ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySig", six * NGENES, NGENES);
        double const * sig_p = sig.data();
        unsigned int const * ranks_p = (unsigned int const *)ranks.data();
        // cache i/o ???

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto & q_up_indices = q_up_indexed[qix];
            auto & q_dn_indices = q_dn_indexed[qix];

            auto const QUSIZE = q_up_indices.size();
            auto const QDSIZE = q_dn_indices.size();

            auto calc_wtks = [&hranks, ranks_p, sig_p](size_type const QSIZE, std::vector<std::uint16_t> const & q_indices)
            {
                double Sum_abs_scores = 0.;

                for (auto & bucket : hranks)
                {
                    bucket.clear();
                }
                for (auto const ix : q_indices)
                {
                    auto abs_score = std::abs(sig_p[ix]);
                    Sum_abs_scores += abs_score;
                    hranks[ranks_p[ix] / RANK_PER_BUCKET].emplace_back((rank_score_t){ix, ranks_p[ix]});
                }
                for (auto & bucket : hranks)
                {
                    std::uint32_t * bucket_p = (std::uint32_t *)bucket.data();
                    std::sort(bucket_p, bucket_p + bucket.size(),
                        [](std::uint32_t p, std::uint32_t q)
                        { return p < q;});
                }

                double const penalty = -1. / (NGENES - QSIZE);
                double wtks = 0.;
                double _acc = 0.;
                double _min = 0.;
                double _max = 0.;
                int prev = 0; // rank indices are 1-based

                for (auto const & bucket : hranks)
                {
                    for (auto const & rs : bucket)
                    {
                        _acc += (rs.rs.rank - prev - 1) * penalty;
                        _min = std::min(_min, _acc);

                        prev = rs.rs.rank;

                        _acc += std::abs(sig_p[rs.rs.score_ix]) / Sum_abs_scores;
                        _max = std::max(_max, _acc);
                    }
                }
                if (_max > std::abs(_min))
                {
                    wtks = _max;
                }
                else
                {
                    wtks = _min;
                }
                return wtks;
            };

            auto const wtks_up = calc_wtks(QUSIZE, q_up_indices);
            auto const wtks_dn = calc_wtks(QDSIZE, q_dn_indices);

            auto wtks = (wtks_dn * wtks_up < 0.) ? (wtks_up - wtks_dn) / 2 : 0.;
            ret.push_back(wtks);

            if (UNLIKELY(m_gt != nullptr))
            {
                auto ref = (*m_gt)[(six - NSKIP) * NQRY + qix];
                if (std::abs(ref - wtks) >= 0.001)
                {
                    std::cout << "Sig: " << six + 1 << ", Query: " << qix + 1 \
                        << " Ref: " << ref << " WTKS: " << wtks << " up/dn " << wtks_up << ' ' << wtks_dn << '\n';
                }
            }
        }
    }

    return ret;
}

#endif /* SRC_CMAP2_HPP_ */
