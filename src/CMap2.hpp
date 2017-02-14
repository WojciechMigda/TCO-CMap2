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

#include "chunked.hpp"
#include "CMAPLib.hpp"
#include "query_parser.hpp"

#include "simd.hpp"

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
#include <climits>

#undef NSIG

enum { PARAM_ROWS_PER_CHUNK_INT = 40000 };
enum { PARAM_ROWS_PER_CHUNK_DBL = 20000 };

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
    {
        m_gene_to_idx.reserve(NGENES);
    }

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

    for (int gix = 0u; gix < genes.size(); ++gix)
    {
        m_gene_to_idx.emplace(genes[gix], gix);
    }

#ifdef CREATE_INDEX_FILES
    auto time0 = timestamp();

    std::cout << "[init] Creating ranksBySigInv chunks" << std::endl;
    create_index_files(PATH + "ranksBySig", "ranksBySigInv", NSIG, NGENES, PARAM_ROWS_PER_CHUNK_INT);
    std::cout << "[init] Elapsed " << timestamp() - time0 << " secs" << std::endl; // 413 sec

    std::cout << "[init] Creating scoresBySigSortedAbsF32 chunks" << std::endl;
    create_score_files(PATH + "scoresBySig", "scoresBySigSortedAbsF32", NSIG, NGENES, PARAM_ROWS_PER_CHUNK_DBL);
    std::cout << "[init] Elapsed " << timestamp() - time0 << " secs" << std::endl; // 1281 sec
#endif

    return 0;
}

template<typename Tp, size_type NCOLS, size_type NROWS, size_type ROW_MAX, size_type ROWS_PER_CHUNK>
struct IOProxy
{
    IOProxy(std::string fn) : m_fn(fn), m_row(INT_MAX), m_v(NCOLS * NROWS)
    {}

    Tp * read_row(size_type row)
    {
        if ((row < m_row) || (row >= (m_row + NROWS)))
        {
            auto start_row = row - row % NROWS;
            auto actual_rows = std::min(ROW_MAX, start_row + NROWS) - start_row;
            auto chunk_num = start_row / ROWS_PER_CHUNK;

            m_v = m_cmap_lib.loadFromIntFile(
                chunk_fname(m_fn, chunk_num),
                ((start_row - chunk_num * ROWS_PER_CHUNK) * NCOLS * sizeof (Tp)) / sizeof (int),
                (actual_rows * NCOLS * sizeof (Tp)) / sizeof (int));
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
    IOProxy<std::uint16_t, NGENES, 10000, 476251, PARAM_ROWS_PER_CHUNK_INT> rank_cache("ranksBySigInv");
    IOProxy<float, NGENES, 10000, 476251, PARAM_ROWS_PER_CHUNK_DBL> score_cache("scoresBySigSortedAbsF32");

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
            std::vector<query_indexed_t> ret;
            ret.reserve(qry.size());

            for (auto const & s : qry)
            {
                auto q_parsed = parse_query(s);
                std::transform(q_parsed.begin(), q_parsed.end(), q_parsed.begin(),
                    [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });
                ret.emplace_back(q_parsed.cbegin(), q_parsed.cend());
            }
            return ret;
        };

    auto const q_up_indexed = q_indexer(q_up);
    auto const q_dn_indexed = q_indexer(q_dn);

//    typedef struct
//    {
//        std::vector<score_index_t> ixs;
////        __v4sf min_max;
//    } query_stream_t;

    using stream_index_t = unsigned int;

    using query_stream_t = std::vector<score_index_t>;

//    std::vector<query_stream_t> q_up_streams(NQRY);
//    std::vector<query_stream_t> q_dn_streams(NQRY);
    // both up and down queries
    std::vector<query_stream_t> q_streams(2 * NQRY);
    std::vector<float> mins(2 * NQRY);
    std::vector<float> maxs(2 * NQRY);

    std::vector<query_stream_t *> gene_buckets[NGENES];
    //std::vector<stream_index_t> gene_buckets[NGENES];

    for (stream_index_t qix = 0u; qix < NQRY; ++qix)
    {
        q_streams[2 * qix + 0]./*ixs.*/reserve(200);
        q_streams[2 * qix + 1]./*ixs.*/reserve(200);

        for (auto const gene_ix : q_up_indexed[qix])
        {
            //gene_buckets[gene_ix].push_back(2 * qix + 0);
            gene_buckets[gene_ix].push_back(&q_streams[2 * qix + 0]);
        }
        for (auto const gene_ix : q_dn_indexed[qix])
        {
            //gene_buckets[gene_ix].push_back(2 * qix + 1);
            gene_buckets[gene_ix].push_back(&q_streams[2 * qix + 1]);
        }
    }

    std::vector<stream_index_t> proc_order(q_streams.size());
    std::iota(proc_order.begin(), proc_order.end(), 0);

    // sort query streams by ascending score indices vector size
    // TODO
    std::sort(proc_order.begin(), proc_order.end(),
        [&q_up_indexed, &q_dn_indexed](stream_index_t p, stream_index_t q)
        {
            std::vector<query_stream_t> const & sp = p % 2 ? q_dn_indexed : q_up_indexed;
            std::vector<query_stream_t> const & sq = q % 2 ? q_dn_indexed : q_up_indexed;

            p /= 2;
            q /= 2;

            return sp[p].size() < sq[q].size();
        });

    for (auto six = NSKIP; six < NSIG; ++six)
    {

        auto inv_ranks = rank_cache.read_row(six);

        auto sigs = score_cache.read_row(six);

        for (auto qix = 0u; qix < 2 * NQRY; ++qix)
        {
            q_streams[qix]./*ixs.*/clear();
        }

        for (auto gix = 0; gix < NGENES; ++gix)
        {
            auto inv_rank = inv_ranks[gix];

            for (auto q_stream_ix : gene_buckets[inv_rank])
            {
                //q_streams[q_stream_ix].ixs.push_back(gix);
                q_stream_ix->/*ixs.*/push_back(gix);
            }
        }

        auto const N4 = NQRY / 2;
        for (auto b4_ix = 0u; b4_ix < N4; ++b4_ix)
        {
//            auto qsix1 = 4 * b4_ix + 0;
//            auto qsix2 = 4 * b4_ix + 1;
//            auto qsix3 = 4 * b4_ix + 2;
//            auto qsix4 = 4 * b4_ix + 3;
            auto qsix1 = proc_order[4 * b4_ix + 0];
            auto qsix2 = proc_order[4 * b4_ix + 1];
            auto qsix3 = proc_order[4 * b4_ix + 2];
            auto qsix4 = proc_order[4 * b4_ix + 3];

            calc_min_max_4(
                q_streams[qsix1], q_streams[qsix2],
                q_streams[qsix3], q_streams[qsix4],
                sigs,
                NGENES,
                mins[qsix1], mins[qsix2], mins[qsix3], mins[qsix4],
                maxs[qsix1], maxs[qsix2], maxs[qsix3], maxs[qsix4]);
        }
        if (UNLIKELY(NQRY % 2))
        {
            auto qsix1 = proc_order[4 * N4 + 0];
            auto qsix2 = proc_order[4 * N4 + 1];

            calc_min_max_2(
                q_streams[qsix1], q_streams[qsix2],
                sigs,
                NGENES,
                mins[qsix1], mins[qsix2],
                maxs[qsix1], maxs[qsix2]);
        }

        for (auto b4_ix = 0u; b4_ix < N4; ++b4_ix)
        {
            __v4sf _min = _mm_load_ps(&mins[4 * b4_ix]);
            __v4sf _max = _mm_load_ps(&maxs[4 * b4_ix]);
            auto vmask = _mm_cmpgt_ps(_max, _mm_and_ps(abs_mask(), _min));
            auto wtks = _mm_or_ps(_mm_and_ps(vmask, _max), _mm_andnot_ps(vmask, _min));

            double _wtks = 0.;
            if (std::signbit(wtks[0]) != std::signbit(wtks[1]))
            {
                _wtks = (wtks[0] - wtks[1]) / 2;
            }
            ret.push_back(_wtks);
                        if (UNLIKELY(m_gt != nullptr))
                        {
                            auto ref = (*m_gt)[(six - NSKIP) * NQRY + b4_ix * 2 + 0];
                            if (std::abs(ref - _wtks) >= 0.001)
                            {
                                std::cout << "Sig: " << six + 1 << ", Query: " << b4_ix * 2 + 0 + 1 \
                                    << " Ref: " << ref << " WTKS: " << _wtks << " up/dn " << wtks[0] << ' ' << wtks[1] << '\n';
                            }
                        }


            _wtks = 0.;
            if (std::signbit(wtks[2]) != std::signbit(wtks[3]))
            {
                _wtks = (wtks[2] - wtks[3]) / 2;
            }
            ret.push_back(_wtks);
                        if (UNLIKELY(m_gt != nullptr))
                        {
                            auto ref = (*m_gt)[(six - NSKIP) * NQRY + b4_ix * 2 + 1];
                            if (std::abs(ref - _wtks) >= 0.001)
                            {
                                std::cout << "Sig: " << six + 1 << ", Query: " << b4_ix * 2 + 0 + 1 \
                                    << " Ref: " << ref << " WTKS: " << _wtks << " up/dn " << wtks[2] << ' ' << wtks[3] << '\n';
                            }
                        }
        }
        // TODO dla nieparzystego
        if (UNLIKELY(NQRY % 2))
        {
            auto _min = (__v4sf)_mm_cvtsi64_si128(pack_2f(mins[4 * N4 + 0], mins[4 * N4 + 1]));
            auto _max = (__v4sf)_mm_cvtsi64_si128(pack_2f(maxs[4 * N4 + 0], maxs[4 * N4 + 1]));

            auto vmask = _mm_cmpgt_ps(_max, _mm_and_ps(abs_mask(), _min));
            auto wtks = _mm_or_ps(_mm_and_ps(vmask, _max), _mm_andnot_ps(vmask, _min));

            double _wtks = 0.;
            if (std::signbit(wtks[0]) != std::signbit(wtks[1]))
            {
                _wtks = (wtks[0] - wtks[1]) / 2;
            }
            ret.push_back(_wtks);
                        if (UNLIKELY(m_gt != nullptr))
                        {
                            auto ref = (*m_gt)[(six - NSKIP) * NQRY + N4 * 2 + 0];
                            if (std::abs(ref - _wtks) >= 0.001)
                            {
                                std::cout << "Sig: " << six + 1 << ", Query: " << N4 * 2 + 0 + 1 \
                                    << " Ref: " << ref << " WTKS: " << _wtks << " up/dn " << wtks[0] << ' ' << wtks[1] << '\n';
                            }
                        }

        }


        for (auto qix = 0u; qix < NQRY; ++qix)
        {

//            double wtks_up;
//            double wtks_dn;
//
//            if (maxs[2 * qix + 0] > std::abs(mins[2 * qix + 0]))
//            {
//                wtks_up = maxs[2 * qix + 0];
//            }
//            else
//            {
//                wtks_up = mins[2 * qix + 0];
//            }
//            if (maxs[2 * qix + 1] > std::abs(mins[2 * qix + 1]))
//            {
//                wtks_dn = maxs[2 * qix + 1];
//            }
//            else
//            {
//                wtks_dn = mins[2 * qix + 1];
//            }
//            auto wtks = (wtks_dn * wtks_up < 0.) ? (wtks_up - wtks_dn) / 2 : 0.;
//            ret.push_back(wtks);
        }

        for (auto qix = 0u; qix < NQRY; ++qix)
        {

//            double wtks_up_dn[2];
//            //calc_wtks_2(q_streams[proc_order[2 * qix]].ixs, q_streams[proc_order[2 * qix + 1]].ixs, sigs, NGENES, wtks_up_dn);
//            calc_wtks_2(q_streams[2 * qix]/*.ixs*/, q_streams[2 * qix + 1]/*.ixs*/, sigs, NGENES, wtks_up_dn);
//            auto wtks_up = wtks_up_dn[0];
//            auto wtks_dn = wtks_up_dn[1];
//
//            auto wtks = (wtks_dn * wtks_up < 0.) ? (wtks_up - wtks_dn) / 2 : 0.;
//            ret.push_back(wtks);
//
//            if (UNLIKELY(m_gt != nullptr))
//            {
//                auto ref = (*m_gt)[(six - NSKIP) * NQRY + qix];
//                if (std::abs(ref - wtks) >= 0.001)
//                {
//                    std::cout << "Sig: " << six + 1 << ", Query: " << qix + 1 \
//                        << " Ref: " << ref << " WTKS: " << wtks << " up/dn " << wtks_up << ' ' << wtks_dn << '\n';
//                }
//            }
        }
    }

    return ret;
}

#endif /* SRC_CMAP2_HPP_ */
