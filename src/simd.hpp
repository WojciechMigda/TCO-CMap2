/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: simd.hpp
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
 * 2017-02-11   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_SIMD_HPP_
#define SRC_SIMD_HPP_

#include "unsafevector.hpp"

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <cmath>

#include <xmmintrin.h>
#include <emmintrin.h>

using query_stream_t = unsafe_vector<std::uint16_t>;

static inline
std::uint64_t pack_2f(float fa, float fb)
{
    typedef union
    {
        float f;
        std::uint32_t i;
    } xf_t;
    xf_t a = {fa};
    xf_t b = {fb};

    return ((std::uint64_t)b.i << 32) | a.i;
}

static inline
__m128 abs_mask(void)
{
    __m128i minus1 = _mm_set1_epi32(-1);
    return _mm_castsi128_ps(_mm_srli_epi32(minus1, 1));
}


// invariant: streams are passed in ascending order of their sizes.
static inline
void calc_min_max_4(
    query_stream_t const & stream1,
    query_stream_t const & stream2,
    query_stream_t const & stream3,
    query_stream_t const & stream4,
    float const * sigs,
    unsigned int NGENES,
    float & omin1,
    float & omin2,
    float & omin3,
    float & omin4,
    float & omax1,
    float & omax2,
    float & omax3,
    float & omax4)
{
    __v4sf vsums = (__v4sf)_mm_setzero_ps();
    auto min_run34 = std::min(stream3.size(), stream4.size());
    auto min_run234 = std::min(stream2.size(), min_run34);
    auto min_run1234 = std::min(stream1.size(), min_run234);
    {
        auto stix = 0u;
        // doing 1,2,3,4
        for (; stix < min_run1234; ++stix)
        {
            auto packed_scores12 = pack_2f(sigs[stream1[stix]], sigs[stream2[stix]]);
            auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
            vsums += (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        }
        // 1 done, doing 2,3,4
        for (; stix < min_run234; ++stix)
        {
            auto packed_scores12 = pack_2f(0, sigs[stream2[stix]]);
            auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
            vsums += (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        }
        // 1,2 done, doing 3,4
        for (; stix < min_run34; ++stix)
        {
            auto packed_scores12 = 0;// pack_2f(0, 0);
            auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
            vsums += (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        }
        // 1,2,3 done, doing 4
        for (; stix < stream4.size(); ++stix)
        {
            auto packed_scores12 = 0;//pack_2f(0, 0);
            auto packed_scores34 = pack_2f(0, sigs[stream4[stix]]);
            vsums += (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        }
    }

    __v4su QSIZE = (__v4su)_mm_set_epi32(stream4.size(), stream3.size(), stream2.size(), stream1.size());
    // TODO, czy tak sie laduje 1-ki ?, ZERO + 1 ?
    __v4sf ONES = (__v4sf){1., 1., 1., 1.};
    __v4sf penalty = _mm_div_ps(ONES, (__v4sf)(_mm_cvtepi32_ps((__m128i)(QSIZE - NGENES))));
    __v4sf divisor = _mm_div_ps(ONES, vsums);

    __v4sf _acc = (__v4sf)_mm_setzero_ps();
    __v4sf _min = (__v4sf)_mm_setzero_ps();
    __v4sf _max = (__v4sf)_mm_setzero_ps();
    __v4su prev = (__v4su)_mm_setzero_si128();

    auto stix = 0u;
    // doing 1,2,3,4
    for (; stix < min_run1234; ++stix)
    {
        auto packed_scores12 = pack_2f(sigs[stream1[stix]], sigs[stream2[stix]]);
        auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
        auto scores = (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        __v4su ix = (__v4su)_mm_set_epi32(stream4[stix], stream3[stix], stream2[stix], stream1[stix]);
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    // 1 done, doing 2,3,4
    penalty[0] = 0.0;
    divisor[0] = 0.0;
    for (; stix < min_run234; ++stix)
    {
        auto packed_scores12 = pack_2f(0, sigs[stream2[stix]]);
        auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
        auto scores = (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        __v4su ix = (__v4su)_mm_set_epi32(stream4[stix], stream3[stix], stream2[stix], 0);
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    // 1,2 done, doing 3,4
    penalty[1] = 0.0;
    divisor[1] = 0.0;
    for (; stix < min_run34; ++stix)
    {
        auto packed_scores12 = 0;//pack_2f(0, 0);
        auto packed_scores34 = pack_2f(sigs[stream3[stix]], sigs[stream4[stix]]);
        auto scores = (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        __v4su ix = (__v4su)_mm_set_epi32(stream4[stix], stream3[stix], 0, 0);
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    // 1,2,3 done, doing 4
    penalty[2] = 0.0;
    divisor[2] = 0.0;
    for (; stix < stream4.size(); ++stix)
    {
        auto packed_scores12 = 0;//pack_2f(0, 0);
        auto packed_scores34 = pack_2f(0, sigs[stream4[stix]]);
        auto scores = (__v4sf)_mm_set_epi64x(packed_scores34, packed_scores12);
        __v4su ix = (__v4su)_mm_set_epi32(stream4[stix], 0, 0, 0);
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }

    union fi
    {
        float f;
        int i;
    };
    _mm_stream_si32((int *)&omin1, (union fi){_min[0]}.i);
    _mm_stream_si32((int *)&omin2, (union fi){_min[1]}.i);
    _mm_stream_si32((int *)&omin3, (union fi){_min[2]}.i);
    _mm_stream_si32((int *)&omin4, (union fi){_min[3]}.i);

    _mm_stream_si32((int *)&omax1, (union fi){_max[0]}.i);
    _mm_stream_si32((int *)&omax2, (union fi){_max[1]}.i);
    _mm_stream_si32((int *)&omax3, (union fi){_max[2]}.i);
    _mm_stream_si32((int *)&omax4, (union fi){_max[3]}.i);

//    omin1 = _min[0];
//    omin2 = _min[1];
//    omin3 = _min[2];
//    omin4 = _min[3];
//    omax1 = _max[0];
//    omax2 = _max[1];
//    omax3 = _max[2];
//    omax4 = _max[3];
}


static inline
void calc_wtks_2(
    std::vector<std::uint16_t> const & up_stream,
    std::vector<std::uint16_t> const & dn_stream,
    float const * sigs,
    unsigned int NGENES,
    double (&o_wtks)[2])
{
    __v4sf vsums = (__v4sf)_mm_setzero_ps();
    {
        auto stix = 0u;
        for (; stix < std::min(up_stream.size(), dn_stream.size()); ++stix)
        {
            auto packed_scores = pack_2f(sigs[up_stream[stix]], sigs[dn_stream[stix]]);
            vsums += (__v4sf)_mm_cvtsi64_si128(packed_scores);
        }
        for (; stix < up_stream.size(); ++stix)
        {
            vsums[0] += sigs[up_stream[stix]];
        }
        for (; stix < dn_stream.size(); ++stix)
        {
            vsums[1] += sigs[dn_stream[stix]];
        }
    }
    vsums += (__v4sf){0., 0., 1., 1.}; // to avoid division by 0

    __v4su QSIZE = (__v4su)_mm_set_epi32(0, 0, dn_stream.size(), up_stream.size());
    // TODO, czy tak sie laduje 1-ki ?, ZERO + 1 ?
    __v4sf ONES = (__v4sf){1., 1., 1., 1.};
    __v4sf penalty = _mm_div_ps(ONES, (__v4sf)(_mm_cvtepi32_ps((__m128i)(QSIZE - NGENES))));
    __v4sf divisor = _mm_div_ps(ONES, vsums);

    __v4sf _acc = (__v4sf)_mm_setzero_ps();
    __v4sf _min = (__v4sf)_mm_setzero_ps();
    __v4sf _max = (__v4sf)_mm_setzero_ps();
    __v4su prev = (__v4su)_mm_setzero_si128();

    auto stix = 0u;
    for (; stix < std::min(up_stream.size(), dn_stream.size()); ++stix)
    {
        auto scores = (__v4sf)_mm_cvtsi64_si128(pack_2f(sigs[up_stream[stix]], sigs[dn_stream[stix]]));
        __v4su ix = (__v4su)_mm_set_epi32(0, 0, dn_stream[stix], up_stream[stix]);
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    for (; stix < up_stream.size(); ++stix)
    {
        auto scores = (__v4sf)_mm_cvtsi64_si128(pack_2f(sigs[up_stream[stix]], 0));
        auto ix = prev;
        ix[0] = up_stream[stix];
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    for (; stix < dn_stream.size(); ++stix)
    {
        auto scores = (__v4sf)_mm_cvtsi64_si128(pack_2f(0, sigs[dn_stream[stix]]));
        auto ix = prev;
        ix[1] = dn_stream[stix];
        _acc += penalty * (__v4sf)_mm_cvtepi32_ps((__m128i)(ix - prev));
        _min = _mm_min_ps(_min, _acc);

        prev = ix + 1;

        _acc += scores * divisor;
        _max = _mm_max_ps(_max, _acc);
    }
    auto vmask = _mm_cmpgt_ps(_max, _mm_and_ps(abs_mask(), _min));
    auto wtks = _mm_or_ps(_mm_and_ps(vmask, _max), _mm_andnot_ps(vmask, _min));
    o_wtks[0] = wtks[0];
    o_wtks[1] = wtks[1];
}


#endif /* SRC_SIMD_HPP_ */
