/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: run_seq.cpp
 *
 * Description:
 *      Single thread with everything sequential, for reference building
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2017-02-13   wm              Initial version
 *
 ******************************************************************************/

#include "query_parser.hpp"
#include "simd.hpp"
#include "likely.h"

#include "cpplinq.hpp"

#include <iostream>
#include <cassert>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <cstdint>
#include <algorithm>
#include <thread>

#include <xmmintrin.h>
#include <emmintrin.h>


enum {NGENES = 10174};
enum {NSIGS = 476251};

enum { PARAM_ROWS_PER_CHUNK_INT = 40000 };
enum { PARAM_ROWS_PER_CHUNK_DBL = 20000 };


using size_type = std::size_t;

using stream_index_t = unsigned int;
using score_index_t = std::uint16_t;
using gene_index_t = std::uint16_t;
using query_indexed_t = std::vector<gene_index_t>;

std::vector<std::string>
read_file_csv(const char * fname)
{
    std::ifstream fcsv(fname);
    std::vector<std::string> vcsv;

    for (std::string line; std::getline(fcsv, line); /* nop */)
    {
        vcsv.push_back(line);
    }
    fcsv.close();

    return vcsv;
}


auto read_genes_to_indices_map = [](char const * fname)
{
    std::unordered_map<int, unsigned int> gene_to_idx;

    FILE * ifile = fopen(fname, "rb");
    if (ifile == nullptr)
    {
        std::cout << "Failed to open " << fname << std::endl;
    }
    else
    {
        fseek(ifile, 0, SEEK_END);
        auto const fsize = ftell(ifile);
        auto const ngenes = fsize / sizeof (int);
        gene_to_idx.reserve(ngenes);

        int * genes = (int *)malloc(fsize);
        assert(genes != nullptr);
        fseek(ifile, 0, SEEK_SET);
        auto const nread = fread(genes, sizeof (int), ngenes, ifile);
        assert(nread == ngenes);
        fclose(ifile);

        std::cout << "Read " << nread << " genes ids" << std::endl;

        for (auto gix = 0u; gix < ngenes; ++gix)
        {
            gene_to_idx.emplace(genes[gix], gix);
        }

        free(genes);
    }

    return gene_to_idx;
};

template<typename Tp>
void load_from_file(
    Tp * obuf_p,
    std::string const & fname,
    size_type pos,
    size_type nelem)
{
    FILE * ifile = fopen(fname.c_str(), "rb");

    if (ifile)
    {
        fseek(ifile, pos * sizeof (Tp), SEEK_SET);
        auto nread = fread(obuf_p, sizeof (Tp), nelem, ifile);
        assert(nread == nelem);
        fclose(ifile);
    }
    else
    {
        std::cout << "Failed to open file " << fname << std::endl;
    }
}

std::string chunk_fname(std::string const & base, std::size_t num)
{
    return base + ".s" + std::to_string(num + 1);
}

template<typename Tp, size_type NCOLS, size_type NROWS, size_type ROW_MAX, size_type ROWS_PER_CHUNK>
struct IOProxy
{
    IOProxy(std::string fn) : m_fn(fn), m_row(INT_MAX)
    {
        m_v = static_cast<Tp *>(malloc(NCOLS * NROWS * sizeof (Tp)));
        assert(m_v != nullptr);
    }
    ~IOProxy()
    {
        free(m_v);
    }

    Tp * read_row(size_type row)
    {
        if ((row < m_row) || (row >= (m_row + NROWS)))
        {
            auto start_row = row - row % NROWS;
            auto actual_rows = std::min(ROW_MAX, start_row + NROWS) - start_row;
            auto chunk_num = start_row / ROWS_PER_CHUNK;

            load_from_file<Tp>(
                m_v,
                chunk_fname(m_fn, chunk_num),
                (start_row - chunk_num * ROWS_PER_CHUNK) * NCOLS,
                actual_rows * NCOLS);

            m_row = start_row;
        }

        return m_v + (row - m_row) * NCOLS;
    }

    std::string m_fn;
    size_type m_row;
    Tp * m_v;
//    std::vector<int> m_v;
//    CMAPLib m_cmap_lib;
};


typedef struct
{
    std::vector<query_indexed_t> const & q_up_indexed;
    std::vector<query_indexed_t> const & q_dn_indexed;

    size_type SIG_BEGIN;
    size_type SIG_END;

    char const * o_wtks_fname;
} worker_ctx_t;

void worker_fn(worker_ctx_t const & ctx)
{
    IOProxy<std::uint16_t, NGENES,
        //10000,
        500,
        NSIGS, PARAM_ROWS_PER_CHUNK_INT> rank_cache("ranksBySigInv");
    IOProxy<float, NGENES,
        //10000,
        500,
        NSIGS, PARAM_ROWS_PER_CHUNK_DBL> score_cache("scoresBySigSortedAbsF32");

    using query_stream_t = std::vector<score_index_t>;

    auto const NQRY = ctx.q_up_indexed.size();
    auto const NQSTREAMS = 2 * NQRY;

    std::vector<query_stream_t> q_streams(NQSTREAMS);

    std::vector<query_stream_t *> gene_buckets[NGENES];

    for (stream_index_t qix = 0u; qix < NQRY; ++qix)
    {
        q_streams[2 * qix + 0].reserve(200);
        q_streams[2 * qix + 1].reserve(200);

        for (auto const gene_ix : ctx.q_up_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_streams[2 * qix + 0]);
        }
        for (auto const gene_ix : ctx.q_dn_indexed[qix])
        {
            gene_buckets[gene_ix].push_back(&q_streams[2 * qix + 1]);
        }
    }

    auto mins = static_cast<float *>(malloc(sizeof (float) * NQSTREAMS));
    auto maxs = static_cast<float *>(malloc(sizeof (float) * NQSTREAMS));

    auto ret_p = static_cast<double *>(malloc(sizeof (double) * NQRY * (ctx.SIG_END - ctx.SIG_BEGIN)));
    size_type ret_ix = 0;

    for (auto six = ctx.SIG_BEGIN; six < ctx.SIG_END; ++six)
    {
        auto inv_ranks = rank_cache.read_row(six);
        auto sigs = score_cache.read_row(six);

        for (auto qix = 0u; qix < NQSTREAMS; ++qix)
        {
            q_streams[qix].clear();
        }

        for (auto gix = 0; gix < NGENES; ++gix)
        {
            auto inv_rank = inv_ranks[gix];

            for (auto q_stream_p : gene_buckets[inv_rank])
            {
                q_stream_p->push_back(gix);
            }
        }

        auto const N4 = NQSTREAMS / 4;

        for (auto b4_ix = 0u; b4_ix < N4; ++b4_ix)
        {
            calc_min_max_4(
                q_streams[4 * b4_ix + 0], q_streams[4 * b4_ix + 1],
                q_streams[4 * b4_ix + 2], q_streams[4 * b4_ix + 3],
                sigs,
                NGENES,
                mins[4 * b4_ix + 0], mins[4 * b4_ix + 1],
                mins[4 * b4_ix + 2], mins[4 * b4_ix + 3],
                maxs[4 * b4_ix + 0], maxs[4 * b4_ix + 1],
                maxs[4 * b4_ix + 2], maxs[4 * b4_ix + 3]);

        }
        if (UNLIKELY(NQRY % 2))
        {
            calc_min_max_2(
                q_streams[4 * N4 + 0], q_streams[4 * N4 + 1],
                sigs,
                NGENES,
                mins[4 * N4 + 0], mins[4 * N4 + 1],
                maxs[4 * N4 + 0], maxs[4 * N4 + 1]);
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
            ret_p[ret_ix++] = _wtks;

            _wtks = 0.;
            if (std::signbit(wtks[2]) != std::signbit(wtks[3]))
            {
                _wtks = (wtks[2] - wtks[3]) / 2;
            }
            ret_p[ret_ix++] = _wtks;
        }
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
            ret_p[ret_ix++] = _wtks;
        }

    }

    FILE * ofile = fopen(ctx.o_wtks_fname, "a+b");
    fseek(ofile, ctx.SIG_BEGIN * NQRY * sizeof (double), SEEK_SET);
    fwrite(ret_p, sizeof (double), ret_ix, ofile);
    fclose(ofile);

    free(ret_p);
    free(mins);
    free(maxs);
}


int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        std::cout << "I need three file names: up and down query CSVs, and an output binary" << std::endl;
        return 1;
    }
    char const * q_up_fname = argv[1];
    char const * q_dn_fname = argv[2];
    char const * o_wtks_fname = argv[3];

    ////////////////////////////////////////////////////////////////////////////

    auto const gene_to_idx = read_genes_to_indices_map("../data/genes.bin");

    auto genes_to_indices = [&gene_to_idx](std::vector<std::string> const & qry) -> std::vector<query_indexed_t>
    {
        std::vector<query_indexed_t> ret;
        ret.reserve(qry.size());

        for (auto const & s : qry)
        {
            auto q_parsed = parse_query(s);
            std::transform(q_parsed.begin(), q_parsed.end(), q_parsed.begin(),
                [&gene_to_idx](int gene_id){ return gene_to_idx.at(gene_id); });
            ret.emplace_back(q_parsed.cbegin(), q_parsed.cend());
        }
        return ret;
    };

    ////////////////////////////////////////////////////////////////////////////

    auto query_scrubber = [](std::string const & s)
        {
            auto sec_comma = s.find(',', s.find(',') + 1);
            return s.substr(sec_comma + 1);
        };
    namespace q = ::cpplinq;
    auto const q_up_vstr = q::from(read_file_csv(q_up_fname))
        >> q::select(query_scrubber)
        >> q::to_vector()
        ;

    auto const q_dn_vstr = q::from(read_file_csv(q_dn_fname))
        >> q::select(query_scrubber)
        >> q::to_vector()
        ;

    assert(q_up_vstr.size() == q_dn_vstr.size());
    auto const NQRY = q_up_vstr.size();

    ////////////////////////////////////////////////////////////////////////////

    auto const q_up_indexed = genes_to_indices(q_up_vstr);
    auto const q_dn_indexed = genes_to_indices(q_dn_vstr);

    ////////////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////////////

    worker_ctx_t ctx[1] =
    {
        {q_up_indexed, q_dn_indexed},
    };
    ctx[0].SIG_BEGIN = 0;
    ctx[0].SIG_END = 476251;
//    ctx[0].SIG_END = 100000;
    ctx[0].o_wtks_fname = o_wtks_fname;

    std::thread worker(worker_fn, ctx[0]);
    worker.join();

    return 0;
}
