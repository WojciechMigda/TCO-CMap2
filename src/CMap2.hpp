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

    using query_stream_t = std::vector<score_index_t>;

    // TODO zmerdzuj w jedno
    std::vector<query_stream_t> q_up_streams(NQRY);
    std::vector<query_stream_t> q_dn_streams(NQRY);

    std::vector<query_stream_t *> gene_buckets[NGENES];

//    for (auto & b : gene_buckets)
    {
//        b.reserve(6);
    }

    for (auto qix = 0u; qix < NQRY; ++qix)
    {
        q_up_streams[qix].reserve(200);
        q_dn_streams[qix].reserve(200);

        for (auto const gene_ix : q_up_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_up_streams[qix]);
        }
        for (auto const gene_ix : q_dn_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_dn_streams[qix]);
        }
    }

    for (auto six = NSKIP; six < NSIG; ++six)
    {

        auto inv_ranks = rank_cache.read_row(six);

        auto sigs = score_cache.read_row(six);

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            q_up_streams[qix].clear();
            q_dn_streams[qix].clear();
        }

        for (auto gix = 0; gix < NGENES; ++gix)
        {
            auto inv_rank = inv_ranks[gix];

            for (auto q_stream_p : gene_buckets[inv_rank])
            {
                q_stream_p->push_back(gix);
            }
        }

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto wtks_calc = [&sigs](query_stream_t const & q_stream, float const sum)
                {
                    auto const QSIZE = q_stream.size();
                    float const penalty = -1. / (NGENES - QSIZE);
                    float const divisor = 1. / sum;

                    double wtks = 0.;
                    float _acc = 0.;
                    float _min = 0.;
                    float _max = 0.;
                    int prev = 0;

                    for (auto const & sr : q_stream)
                    {
                        auto ix = sr;
                        auto sc = sigs[sr];
                        _acc += (ix - prev) * penalty;
//                        _acc += (ix - prev) / inv_penalty;
                        _min = std::min(_min, _acc);

                        prev = ix + 1;

//                        _acc += sc / sum;
                        _acc += sc * divisor;
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

            float up_sum = 0.;
            float dn_sum = 0.;
            {
                auto stix = 0u;
                for (; stix < std::min(q_up_streams[qix].size(), q_dn_streams[qix].size()); ++stix)
                {
                    up_sum += sigs[q_up_streams[qix][stix]];
                    dn_sum += sigs[q_dn_streams[qix][stix]];
                }
                for (; stix < q_up_streams[qix].size(); ++stix)
                {
                    up_sum += sigs[q_up_streams[qix][stix]];
                }
                for (; stix < q_dn_streams[qix].size(); ++stix)
                {
                    dn_sum += sigs[q_dn_streams[qix][stix]];
                }
            }


            double wtks_up = wtks_calc(q_up_streams[qix], up_sum);
            double wtks_dn = wtks_calc(q_dn_streams[qix], dn_sum);

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
