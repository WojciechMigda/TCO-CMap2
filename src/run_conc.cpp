/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: run_conc.cpp
 *
 * Description:
 *      Fully threaded execution
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
#include "unsafevector.hpp"
#include "memalign.hpp"

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
#include <cstddef>
#include <chrono>
#include <cstring>

#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>

#include <xmmintrin.h>
#include <emmintrin.h>

//#define USE_MMAP
#define USE_PACKED_RANKS

enum {NGENES = 10174};
enum {NSIGS = 476251};
enum { NPKRANKS = 8903 }; // ceil(NGENES * 7 / 8)

enum { PARAM_ROWS_PER_CHUNK_INT = 40000 };
enum { PARAM_ROWS_PER_CHUNK_DBL = 20000 };

enum { XMM_ALIGN = 16 };

using size_type = std::size_t;
using ssize_type = ssize_t;

using stream_index_t = unsigned int;
using score_index_t = std::uint16_t;
using gene_index_t = std::uint16_t;
using query_indexed_t = std::vector<gene_index_t>;

#include <sys/time.h>
timeval timestamp()
{
    timeval tv;
    gettimeofday(&tv, NULL);

    return tv;
}

int set_affinity(std::thread::native_handle_type tid, unsigned int cpuid)
{
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpuid, &cpuset);
    int rc = pthread_setaffinity_np(tid, sizeof (cpu_set_t), &cpuset);
    if (rc != 0)
    {
        std::cout << "pthread_setaffinity_np failed for " << cpuid << std::endl;
    }
    return rc;
}

std::vector<std::string>
read_file_csv(const char * fname)
{
    std::ifstream fcsv(fname);
    if (!fcsv.is_open())
    {
        std::cout << "Failed to open CSV query file " << fname << std::endl;
        exit(1);
    }

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

    if (UNLIKELY(ifile == nullptr))
    {
        std::cout << "Failed to open genes file " << fname << std::endl;
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
    int idesc,
    size_type pos,
    size_type nelem)
{
    //auto t0 = std::chrono::high_resolution_clock::now();
    //fseek(ifile, pos * sizeof (Tp), SEEK_SET);
    char * op = (char *)obuf_p;
    auto to_read = sizeof (Tp) * nelem;
    ssize_type nread = 0;
    do
    {
        nread = read(idesc, op, to_read);
        if (UNLIKELY(nread == -1))
        {
            assert(errno == EINTR);
            continue;
        }
        else
        {
            to_read -= nread;
            op += nread;
        }
    } while (to_read != 0);
    //std::cout << "Read2 t= " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << std::endl;
    assert(nread == sizeof (Tp) * nelem);
}


template <typename Tp>
void wtks_saver(int odesc, Tp const * ibuf_p, size_type pos, size_type n)
{
    char const * ip = (char const *)ibuf_p;
    auto to_write = sizeof (Tp) * n;
    ssize_type nwr = 0;
    do
    {
        nwr = write(odesc, ip, to_write);
        if (UNLIKELY(nwr == -1))
        {
            assert(errno == EINTR);
            continue;
        }
        else
        {
            to_write -= nwr;
            ip += nwr;
        }
    } while (to_write != 0);

    //std::cout << "Wrote " << n << " records at " << pos << std::endl;
}


using query_stream_t = unsafe_vector<score_index_t>;

typedef struct
{
    std::vector<query_indexed_t> const * q_up_indexed_p;
    std::vector<query_indexed_t> const * q_dn_indexed_p;

    unsigned int id;

    size_type SIG_BEGIN;
    size_type SIG_END;

    float * sigs;
    std::uint16_t * ranks;
    float * mins;
    float * maxs;
    float * owtks;
    std::vector<query_stream_t> q_streams;
    std::vector<stream_index_t> proc_order;
    std::vector<query_stream_t *> gene_buckets[NGENES];
} worker_ctx_t;

void worker(worker_ctx_t * ctx_p)
{
    auto rank_unpk = [](std::uint16_t const * p, std::size_t ix) -> std::uint16_t
        {
            if (LIKELY(ix % 8 != 7))
            {
                return p[ix - ix / 8] & 0x3FFF;
            }
            else
            {
                __v8hu v8 = (__v8hu)_mm_loadu_si128((__m128i *)&p[ix - ix / 8 - 7]);
                std::uint16_t hi = _mm_movemask_epi8((__m128i)v8);
                hi &= 0x2AAA;
                std::uint16_t lo = _mm_movemask_epi8(_mm_slli_epi16((__m128i)v8, 1));
                lo &= 0x2AAA;
                return hi | (lo >> 1);
            }
        };

    auto & ctx = *ctx_p;
    //std::cout << "Worker " << ctx.id << ", sigs=[" << ctx.SIG_BEGIN << ',' << ctx.SIG_END << ')' << std::endl;
    auto & q_up_indexed = *ctx.q_up_indexed_p;

    auto const NQRY = q_up_indexed.size();
    auto const NQSTREAMS = 2 * NQRY;

    auto & q_streams = ctx.q_streams;
    auto const & gene_buckets = ctx.gene_buckets;
    auto const & proc_order = ctx.proc_order;

    auto mins = ctx.mins;
    auto maxs = ctx.maxs;

    size_type ret_ix = 0;

    for (auto six = ctx.SIG_BEGIN; six < ctx.SIG_END; ++six)
    {
#ifdef USE_PACKED_RANKS
        auto inv_ranks = &ctx.ranks[(six - ctx.SIG_BEGIN) * NPKRANKS];
#else
        auto inv_ranks = &ctx.ranks[(six - ctx.SIG_BEGIN) * NGENES];
#endif
        auto sigs = &ctx.sigs[(six - ctx.SIG_BEGIN) * NGENES];

        for (auto qix = 0u; qix < NQSTREAMS; ++qix)
        {
            q_streams[qix].clear();
        }

        for (auto gix = 0; gix < NGENES; ++gix)
        {
#ifdef USE_PACKED_RANKS
            auto inv_rank = rank_unpk(inv_ranks, gix);
#else
            auto inv_rank = inv_ranks[gix];
#endif

            for (auto q_stream_p : gene_buckets[inv_rank])
            {
                q_stream_p->push_back(gix);
            }
        }

        auto const N4 = NQSTREAMS / 4;

        for (auto b4_ix = 0u; b4_ix < N4; ++b4_ix)
        {
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
            // duplicates, so I can reuse the 4-way version
            auto qsix3 = proc_order[4 * N4 + 0];
            auto qsix4 = proc_order[4 * N4 + 1];

            calc_min_max_4(
                q_streams[qsix1], q_streams[qsix2],
                q_streams[qsix3], q_streams[qsix4],
                sigs,
                NGENES,
                mins[qsix1], mins[qsix2], mins[qsix3], mins[qsix4],
                maxs[qsix1], maxs[qsix2], maxs[qsix3], maxs[qsix4]);
        }


        for (auto b4_ix = 0u; b4_ix < N4; ++b4_ix)
        {
            __v4sf _min = _mm_load_ps(&mins[4 * b4_ix]);
            __v4sf _max = _mm_load_ps(&maxs[4 * b4_ix]);
            auto vmask = _mm_cmpgt_ps(_max, _mm_and_ps(abs_mask(), _min));
            auto wtks = _mm_blendv_ps(_min, _max, vmask);

            float _wtks = 0.;
            auto const signums = _mm_movemask_ps(wtks);
            if (((signums & 3) == 1) || ((signums & 3) == 2))
            {
                _wtks = (wtks[0] - wtks[1]) / 2;
            }
            //ctx.owtks[ret_ix++] = _wtks;
            mm_stream_ss(&ctx.owtks[ret_ix], _wtks);
            ++ret_ix;

            _wtks = 0.;
            if (((signums & 0xC) == 8) || ((signums & 0xC) == 4))
            {
                _wtks = (wtks[2] - wtks[3]) / 2;
            }
            mm_stream_ss(&ctx.owtks[ret_ix], _wtks);
            ++ret_ix;

        }
        if (UNLIKELY(NQRY % 2))
        {
            auto _min = (__v4sf)_mm_cvtsi64_si128(pack_2f(mins[4 * N4 + 0], mins[4 * N4 + 1]));
            auto _max = (__v4sf)_mm_cvtsi64_si128(pack_2f(maxs[4 * N4 + 0], maxs[4 * N4 + 1]));

            auto vmask = _mm_cmpgt_ps(_max, _mm_and_ps(abs_mask(), _min));
            auto wtks = _mm_blendv_ps(_min, _max, vmask);

            float _wtks = 0.;
            auto const signums = _mm_movemask_ps(wtks);
            if (((signums & 3) == 1) || ((signums & 3) == 2))
            {
                _wtks = (wtks[0] - wtks[1]) / 2;
            }
            mm_stream_ss(&ctx.owtks[ret_ix], _wtks);
            ++ret_ix;
        }

    }
//    _mm_sfence();
}

typedef struct
{
    unsigned int cpuid;
    std::vector<query_indexed_t> const * q_up_indexed;
    std::vector<query_indexed_t> const * q_dn_indexed;

    size_type SIG_BEGIN;
    size_type SIG_END;

    char const * o_wtks_fname;
} producer_ctx_t;


void producer(producer_ctx_t const & ctx)
{
    std::cout << "Producer " << ctx.cpuid + 1 << ", sigs=[" << ctx.SIG_BEGIN << ',' << ctx.SIG_END << ')' << std::endl;

    char const * sigs_fn = "../data/scoresBySigSortedAbsF32";
#ifdef USE_PACKED_RANKS
    char const * ranks_fn = "../data/ranksBySigInvPacked";
#else
    char const * ranks_fn = "../data/ranksBySigInv";
#endif

    typedef struct
    {
        std::thread worker;
        worker_ctx_t ctx;
    } job_t;

    job_t jobs[2];
    memset(jobs, 0, sizeof (jobs));
    job_t * jobs_p[2] = {&jobs[0], &jobs[1]};

    enum { BATCH_SZ = 512 };
    auto const NQRY = ctx.q_up_indexed->size();
    auto const NQSTREAMS = 2 * NQRY;


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////

    auto worker_ctx_initialize = [NQRY, NQSTREAMS, &ctx](worker_ctx_t & wctx)
    {
        wctx.sigs = static_cast<float *>(aligned_malloc(XMM_ALIGN, BATCH_SZ * NGENES * sizeof (float)));
#ifdef USE_PACKED_RANKS
        wctx.ranks = static_cast<std::uint16_t *>(aligned_malloc(XMM_ALIGN, BATCH_SZ * NPKRANKS * sizeof (std::uint16_t)));
#else
        wctx.ranks = static_cast<std::uint16_t *>(aligned_malloc(XMM_ALIGN, BATCH_SZ * NGENES * sizeof (std::uint16_t)));
#endif
        wctx.owtks = static_cast<float *>(aligned_malloc(XMM_ALIGN, BATCH_SZ * NQRY * sizeof (float)));

        wctx.q_up_indexed_p = ctx.q_up_indexed;
        wctx.q_dn_indexed_p = ctx.q_dn_indexed;

        wctx.q_streams = std::vector<query_stream_t>(NQSTREAMS);
        for (auto & qs : wctx.q_streams) { qs.reserve(200); }

        for (stream_index_t qix = 0u; qix < NQRY; ++qix)
        {
            for (auto const gene_ix : (*ctx.q_up_indexed)[qix])
            {
                wctx.gene_buckets[gene_ix].push_back(&wctx.q_streams[2 * qix + 0]);
            }
            for (auto const gene_ix : (*ctx.q_dn_indexed)[qix])
            {
                wctx.gene_buckets[gene_ix].push_back(&wctx.q_streams[2 * qix + 1]);
            }
        }

        wctx.proc_order = std::vector<stream_index_t>(wctx.q_streams.size());
        std::iota(wctx.proc_order.begin(), wctx.proc_order.end(), 0);

        // sort query streams by ascending score indices vector size
        std::sort(wctx.proc_order.begin(), wctx.proc_order.end(),
            [&ctx](stream_index_t p, stream_index_t q)
            {
                std::vector<query_indexed_t> const & sp = p % 2 ? (*ctx.q_dn_indexed) : (*ctx.q_up_indexed);
                std::vector<query_indexed_t> const & sq = q % 2 ? (*ctx.q_dn_indexed) : (*ctx.q_up_indexed);

                p /= 2;
                q /= 2;

                return sp[p].size() < sq[q].size();
            });

        wctx.mins = static_cast<float *>(aligned_malloc(XMM_ALIGN, sizeof (float) * NQSTREAMS));
        wctx.maxs = static_cast<float *>(aligned_malloc(XMM_ALIGN, sizeof (float) * NQSTREAMS));
    };

    worker_ctx_initialize(jobs[0].ctx);
    worker_ctx_initialize(jobs[1].ctx);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    // stub worker so I can join it later in the loop
    jobs[1].worker = std::thread([](){});

    auto isig_file = open(sigs_fn, O_RDONLY);
    if (UNLIKELY(isig_file == -1))
    {
        std::cout << "Failed to open signature file " << sigs_fn << std::endl;
    }
    lseek(isig_file, ctx.SIG_BEGIN * NGENES * sizeof (float), SEEK_SET);
    posix_fadvise(isig_file, ctx.SIG_BEGIN * NGENES * sizeof (float), ctx.SIG_END * NGENES * sizeof (float), POSIX_FADV_SEQUENTIAL);
    posix_fadvise(isig_file, ctx.SIG_BEGIN * NGENES * sizeof (float), ctx.SIG_END * NGENES * sizeof (float), POSIX_FADV_NOREUSE);


    auto irank_file = open(ranks_fn, O_RDONLY);
    if (UNLIKELY(irank_file == -1))
    {
        std::cout << "Failed to open ranks file " << irank_file << std::endl;
    }
#ifdef USE_PACKED_RANKS
    lseek(irank_file, ctx.SIG_BEGIN * NPKRANKS * sizeof (std::uint16_t), SEEK_SET);
    posix_fadvise(irank_file, ctx.SIG_BEGIN * NPKRANKS * sizeof (std::uint16_t), ctx.SIG_END * NPKRANKS * sizeof (std::uint16_t), POSIX_FADV_SEQUENTIAL);
    posix_fadvise(irank_file, ctx.SIG_BEGIN * NPKRANKS * sizeof (std::uint16_t), ctx.SIG_END * NPKRANKS * sizeof (std::uint16_t), POSIX_FADV_NOREUSE);
#else
    lseek(irank_file, ctx.SIG_BEGIN * NGENES * sizeof (std::uint16_t), SEEK_SET);
    posix_fadvise(irank_file, ctx.SIG_BEGIN * NGENES * sizeof (std::uint16_t), ctx.SIG_END * NGENES * sizeof (std::uint16_t), POSIX_FADV_SEQUENTIAL);
    posix_fadvise(irank_file, ctx.SIG_BEGIN * NGENES * sizeof (std::uint16_t), ctx.SIG_END * NGENES * sizeof (std::uint16_t), POSIX_FADV_NOREUSE);
#endif

    auto ofile = open(ctx.o_wtks_fname, O_WRONLY);
    if (UNLIKELY(ofile == -1))
    {
        std::cout << "Failed to open output WTKS file " << ctx.o_wtks_fname << std::endl;
    }
    lseek(ofile, ctx.SIG_BEGIN * NQRY * sizeof (float), SEEK_SET);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    for (auto six = ctx.SIG_BEGIN; six < ctx.SIG_END; six += BATCH_SZ)
    {
        auto const nrows = std::min<size_type>(ctx.SIG_END - six, BATCH_SZ);

        load_from_file<float>(
            jobs_p[0]->ctx.sigs, isig_file, six * NGENES, nrows * NGENES);
#ifdef USE_PACKED_RANKS
        load_from_file<std::uint16_t>(
            jobs_p[0]->ctx.ranks, irank_file, six * NPKRANKS, nrows * NPKRANKS);
#else
        load_from_file<std::uint16_t>(
            jobs_p[0]->ctx.ranks, irank_file, six * NGENES, nrows * NGENES);
#endif

        jobs_p[0]->ctx.SIG_BEGIN = six;
        jobs_p[0]->ctx.SIG_END = six + nrows;

        jobs_p[0]->ctx.id = (ctx.cpuid + 1) * 10000 + 1 + (six - ctx.SIG_BEGIN) / BATCH_SZ;

        jobs_p[1]->worker.join();

        jobs_p[0]->worker = std::thread(worker, &jobs_p[0]->ctx);
        set_affinity(jobs_p[0]->worker.native_handle(), ctx.cpuid);

        wtks_saver<float>(ofile, jobs_p[1]->ctx.owtks, jobs_p[1]->ctx.SIG_BEGIN * NQRY, (jobs_p[1]->ctx.SIG_END - jobs_p[1]->ctx.SIG_BEGIN) * NQRY);

        std::swap(jobs_p[0], jobs_p[1]);
    }

    jobs_p[1]->worker.join();
    wtks_saver<float>(ofile, jobs_p[1]->ctx.owtks, jobs_p[1]->ctx.SIG_BEGIN * NQRY, (jobs_p[1]->ctx.SIG_END - jobs_p[1]->ctx.SIG_BEGIN) * NQRY);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////


    close(isig_file);
    close(irank_file);
    close(ofile);

    auto worker_ctx_destroy = [](worker_ctx_t const & wctx)
    {
        free(wctx.sigs);
        free(wctx.ranks);
        free(wctx.owtks);
        free(wctx.mins);
        free(wctx.maxs);
    };

    worker_ctx_destroy(jobs[0].ctx);
    worker_ctx_destroy(jobs[1].ctx);
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

    ////////////////////////////////////////////////////////////////////////////

    auto const q_up_indexed = genes_to_indices(q_up_vstr);
    auto const q_dn_indexed = genes_to_indices(q_dn_vstr);

    ////////////////////////////////////////////////////////////////////////////

    FILE * ofile = fopen(o_wtks_fname, "wb");
    if (ofile == nullptr)
    {
        std::cout << "Cannot create/open for writing: " << o_wtks_fname << std::endl;
        return 1;
    }
    // truncate
    {
        auto len = NSIGS * q_up_indexed.size() * sizeof (float);
        auto out = fileno(ofile);
        if ((fallocate(out, 0, 0, len) == -1) && (errno == EOPNOTSUPP))
        {
            auto ret = ftruncate(out, len);
            (void)ret;
        }
    }

    fclose(ofile);

    ////////////////////////////////////////////////////////////////////////////

    std::cout << "Hardware concurrency " << std::thread::hardware_concurrency() << std::endl;

    auto const NTHR = std::thread::hardware_concurrency();

    std::vector<size_type> ranges(NTHR + 1);
    ranges.front() = 0;
    ranges.back() = NSIGS;
    for (auto rix = 1u; rix < NTHR; ++rix)
    {
        ranges[rix] = NSIGS * rix / NTHR;
    }

    std::vector<producer_ctx_t> ctx(NTHR);
    for (auto cix = 0u; cix < NTHR; ++cix)
    {
        ctx[cix].cpuid = cix;
        ctx[cix].o_wtks_fname = o_wtks_fname;
        ctx[cix].q_up_indexed = &q_up_indexed;
        ctx[cix].q_dn_indexed = &q_dn_indexed;
        ctx[cix].SIG_BEGIN = ranges[cix];
        ctx[cix].SIG_END = ranges[cix + 1];
    }

    std::thread threads[NTHR];
    // TODO testing
    for (auto tix = 0u; tix < NTHR; ++tix)
    //for (auto tix = 1u; tix < 2; ++tix)
    {
        threads[tix] = std::thread(producer, ctx[tix]);

        // set afinity
        set_affinity(threads[tix].native_handle(), tix);
    }

    for (auto tix = 0u; tix < NTHR; ++tix)
    //for (auto tix = 1u; tix < 2; ++tix)
    {
        threads[tix].join();
    }

    return 0;
}
