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

inline
void calc_wtks_2(
    std::vector<std::uint16_t> const & up_stream,
    std::vector<std::uint16_t> const & dn_stream,
    float const * sigs,
    std::size_t NGENES,
    double (&o_wtks)[2])
{
    float up_sum = 0.;
    float dn_sum = 0.;
    {
        auto stix = 0u;
        for (; stix < std::min(up_stream.size(), dn_stream.size()); ++stix)
        {
            up_sum += sigs[up_stream[stix]];
            dn_sum += sigs[dn_stream[stix]];
        }
        for (; stix < up_stream.size(); ++stix)
        {
            up_sum += sigs[up_stream[stix]];
        }
        for (; stix < dn_stream.size(); ++stix)
        {
            dn_sum += sigs[dn_stream[stix]];
        }
    }

    auto const QUSIZE = up_stream.size();
    auto const QDSIZE = dn_stream.size();

    float const up_penalty = -1. / (NGENES - QUSIZE);
    float const dn_penalty = -1. / (NGENES - QDSIZE);

    float const up_divisor = 1. / up_sum;
    float const dn_divisor = 1. / dn_sum;

    double wtks_up = 0.;
    float _acc_up = 0.;
    float _min_up = 0.;
    float _max_up = 0.;
    int prev_up = 0;

    double wtks_dn = 0.;
    float _acc_dn = 0.;
    float _min_dn = 0.;
    float _max_dn = 0.;
    int prev_dn = 0;

    for (auto const & sr : up_stream)
    {
        auto ix = sr;
        auto sc = sigs[sr];
        _acc_up += (ix - prev_up) * up_penalty;
//                        _acc += (ix - prev) / inv_penalty;
        _min_up = std::min(_min_up, _acc_up);

        prev_up = ix + 1;

//                        _acc += sc / sum;
        _acc_up += sc * up_divisor;
        _max_up = std::max(_max_up, _acc_up);
    }
    if (_max_up > std::abs(_min_up))
    {
        wtks_up = _max_up;
    }
    else
    {
        wtks_up = _min_up;
    }

    for (auto const & sr : dn_stream)
    {
        auto ix = sr;
        auto sc = sigs[sr];
        _acc_dn += (ix - prev_dn) * dn_penalty;
//                        _acc += (ix - prev) / inv_penalty;
        _min_dn = std::min(_min_dn, _acc_dn);

        prev_dn = ix + 1;

//                        _acc += sc / sum;
        _acc_dn += sc * dn_divisor;
        _max_dn = std::max(_max_dn, _acc_dn);
    }
    if (_max_dn > std::abs(_min_dn))
    {
        wtks_dn = _max_dn;
    }
    else
    {
        wtks_dn = _min_dn;
    }

    o_wtks[0] = wtks_up;
    o_wtks[1] = wtks_dn;
}


#endif /* SRC_SIMD_HPP_ */
