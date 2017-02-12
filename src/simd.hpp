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

#include <cstdint>
#include <cstddef>
#include <algorithm>
#include <cmath>

#include <xmmintrin.h>
#include <emmintrin.h>

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

static inline
void calc_wtks_2(
    std::vector<std::uint16_t> const & up_stream,
    std::vector<std::uint16_t> const & dn_stream,
    float const * sigs,
    unsigned int NGENES,
    double (&o_wtks)[2])
{
    __v4sf vsums = (__v4sf)_mm_setzero_ps();
//    float sums[4];
    {
        auto stix = 0u;
        for (; stix < std::min(up_stream.size(), dn_stream.size()); ++stix)
        {
            auto packed_scores = pack_2f(sigs[up_stream[stix]], sigs[dn_stream[stix]]);
            vsums += (__v4sf)_mm_cvtsi64_si128(packed_scores);
        }
//        _mm_store_ps(sums, vsums);
        for (; stix < up_stream.size(); ++stix)
        {
//            sums[0] += sigs[up_stream[stix]];
            vsums[0] += sigs[up_stream[stix]];
        }
        for (; stix < dn_stream.size(); ++stix)
        {
//            sums[1] += sigs[dn_stream[stix]];
            vsums[1] += sigs[dn_stream[stix]];
        }
    }
    vsums += (__v4sf){0., 0., 1., 1.}; // to avoid division by 0

    __v4su QSIZE = (__v4su)_mm_set_epi32(0, 0, dn_stream.size(), up_stream.size());
    // TODO, czy tak sie laduje 1-ki ?, ZERO + 1 ?
    __v4sf ONES = (__v4sf){1., 1., 1., 1.};
    __v4sf penalty = _mm_div_ps(ONES, (__v4sf)(_mm_cvtepi32_ps((__m128i)(QSIZE - /*_mm_set1_epi32*/(NGENES)))));
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
    return;

//
//    auto const QUSIZE = up_stream.size();
//    auto const QDSIZE = dn_stream.size();
//
//    float const up_penalty = -1. / (NGENES - QUSIZE);
//    float const dn_penalty = -1. / (NGENES - QDSIZE);
//
//    float const up_divisor = 1. / sums[0];
//    float const dn_divisor = 1. / sums[1];
//
//    double wtks_up = 0.;
//    float _acc_up = 0.;
//    float _min_up = 0.;
//    float _max_up = 0.;
//    int prev_up = 0;
//
//    double wtks_dn = 0.;
//    float _acc_dn = 0.;
//    float _min_dn = 0.;
//    float _max_dn = 0.;
//    int prev_dn = 0;
//
//    for (auto const & sr : up_stream)
//    {
//        auto ix = sr;
//        auto sc = sigs[sr];
//        _acc_up += (ix - prev_up) * up_penalty;
////                        _acc += (ix - prev) / inv_penalty;
//        _min_up = std::min(_min_up, _acc_up);
//
//        prev_up = ix + 1;
//
////                        _acc += sc / sum;
//        _acc_up += sc * up_divisor;
//        _max_up = std::max(_max_up, _acc_up);
//    }
//    if (_max_up > std::abs(_min_up))
//    {
//        wtks_up = _max_up;
//    }
//    else
//    {
//        wtks_up = _min_up;
//    }
//
//    for (auto const & sr : dn_stream)
//    {
//        auto ix = sr;
//        auto sc = sigs[sr];
//        _acc_dn += (ix - prev_dn) * dn_penalty;
////                        _acc += (ix - prev) / inv_penalty;
//        _min_dn = std::min(_min_dn, _acc_dn);
//
//        prev_dn = ix + 1;
//
////                        _acc += sc / sum;
//        _acc_dn += sc * dn_divisor;
//        _max_dn = std::max(_max_dn, _acc_dn);
//    }
//    if (_max_dn > std::abs(_min_dn))
//    {
//        wtks_dn = _max_dn;
//    }
//    else
//    {
//        wtks_dn = _min_dn;
//    }
//
//    o_wtks[0] = wtks_up;
//    o_wtks[1] = wtks_dn;
}


#endif /* SRC_SIMD_HPP_ */
