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
#include <cstring>

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

template<typename Tp, size_type NCOLS, size_type NROWS, size_type ROW_MAX>
struct IOProxy
{
    IOProxy(std::string fn) : m_fn(fn), m_row(INT_MAX), m_v(NCOLS * NROWS)
    {}

    Tp * read_row(size_type row)
    {
        if ((row < m_row) || (row >= (m_row + NROWS)))
        {
            auto start_row = row - row % NROWS;
            auto end_row = std::min(ROW_MAX, start_row + NROWS);

            m_v = m_cmap_lib.loadFromIntFile(
                m_fn,
                (start_row * NCOLS * sizeof (Tp)) / sizeof (int),
                ((end_row - start_row) * NCOLS * sizeof (Tp)) / sizeof (int));
            m_row = start_row;
        }

        return (Tp *)m_v.data() + (row - m_row) * NCOLS;
    }

    std::string m_fn;
    size_type m_row;
    std::vector<int> m_v;
    CMAPLib m_cmap_lib;
};

std::vector<double>
CMAP2Updated::getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn)
{
    IOProxy<std::uint16_t, 10174, 10000, 476251> rank_cache(PATH + "ranksBySigInv");
    IOProxy<float, 10174, 10000, 476251> score_cache(PATH + "scoresBySigSortedAbsF32");

    std::cout << "Will process " << NSIG - NSKIP << " signatures, [" << NSKIP << ',' << NSIG << ")\n";

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

//    typedef union
//    {
//        double f64;
//        std::uint64_t i64;
//    } packed_score_ix_t;
    typedef struct
    {
        float score;
        std::uint16_t ix;
    } packed_score_ix_t;
//    using packed_score_ix_t = std::uint64_t;

    //using score_ix_t = std::pair<double, score_index_t>;
    using score_ix_t = packed_score_ix_t;

    typedef struct
    {
        std::vector<score_ix_t> scores;
        double sum;
    } query_state_t;

    std::vector<query_state_t> q_up_states(NQRY);
    std::vector<query_state_t> q_dn_states(NQRY);

    std::vector<query_state_t *> gene_buckets[NGENES];

    for (auto & b : gene_buckets)
    {
//        b.reserve(6);
    }

    for (auto qix = 0u; qix < NQRY; ++qix)
    {
        q_up_states[qix].scores.reserve(200);
        q_dn_states[qix].scores.reserve(200);

        for (auto const gene_ix : q_up_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_up_states[qix]);
        }
        for (auto const gene_ix : q_dn_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_dn_states[qix]);
        }
    }

    for (auto six = NSKIP; six < NSIG; ++six)
    {

//        auto _ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySigInv", six * NGENES / 2, NGENES / 2);
//        auto inv_ranks = (std::uint16_t const *)_ranks.data();
        auto inv_ranks = rank_cache.read_row(six);
        //assert(memcmp(foo, inv_ranks, NGENES * 2) == 0);

//        auto _sigs = m_cmap_lib.loadFromIntFile(PATH + "scoresBySigSortedAbsF32", six * NGENES, NGENES);
//        auto sigs = (float const *)_sigs.data();
        auto sigs = score_cache.read_row(six);

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
            auto score = sigs[gix];

            for (auto q_state_p : gene_buckets[inv_rank])
            {
                packed_score_ix_t sx;
                sx.score = score;
                sx.ix = gix;
                q_state_p->scores.push_back(sx);
                q_state_p->sum += score;
            }
        }

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto wtks_calc = [](query_state_t const & q_state)
                {
                    auto const QSIZE = q_state.scores.size();
                    float const penalty = -1. / (NGENES - QSIZE);
                    double const inv_penalty = NGENES - QSIZE;
                    double const divisor = 1. / q_state.sum;
                    float const hit_fac = q_state.sum;

                    // less cache misses but more insns..
//                    float const hit_fac = ({
//                        double sum = 0.;
//                        for (auto const & sr : q_state.scores)
//                        {
//                            sum += sr.score;
//                        }
//                        sum;
//                    });

                    double wtks = 0.;
                    float _acc = 0.;
                    float _min = 0.;
                    float _max = 0.;
                    int prev = 0; // rank indices are 1-based

                    for (auto const & sr : q_state.scores)
                    {
                        auto ix = sr.ix;
                        auto sc = sr.score;
                        _acc += (ix - prev - 0) * penalty;
//                        _acc += (ix - prev - 0) / inv_penalty;
                        _min = std::min(_min, _acc);

                        prev = ix + 1;

//                        _acc += sc / q_state.sum;
                        _acc += sc / hit_fac;
//                        _acc += sc * divisor;
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

#endif /* SRC_CMAP2_HPP_ */
