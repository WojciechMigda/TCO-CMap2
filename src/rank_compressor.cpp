/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: rank_compressor.cpp
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
 * 2017-02-15   wm              Initial version
 *
 ******************************************************************************/

#include "likely.h"

#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <cstddef>

#include <emmintrin.h>

#undef NDEBUG

enum { NGENES = 10174 };
enum { NSIGS = 476251 };

int main(int argc, char ** argv)
{
    auto irankf = fopen("ranksBySigInv", "rb");
    if (irankf == nullptr)
    {
        std::cout << "Could not open for reading ranksBySigInv" << std::endl;
        return 1;
    }

    auto orankf = fopen("ranksBySigInvPacked", "wb");
    if (orankf == nullptr)
    {
        std::cout << "Could not open for writing ranksBySigInvPacked" << std::endl;
        return 1;
    }

    std::cout << "Assuming " << NSIGS << " signatures and " << NGENES << " genes" << std::endl;


    ////////////////////////////////////////////////////////////////////////////

    {
        std::uint16_t irank_row[NGENES];
        //enum { NPACKED = 3 };
        //std::uint16_t orank_row[NGENES];
        std::vector<std::uint16_t> orank_row;

        for (auto six = 0u; six < NSIGS; ++six)
        {
            auto nread = fread(irank_row, sizeof (irank_row[0]), NGENES, irankf);
            assert(nread == NGENES);

            for (auto cix = 0u; cix < NGENES; ++cix)
            {
                if (cix % 8 == 7)
                {
                    auto const pix = orank_row.size();
                    orank_row[pix - 1] |= (irank_row[cix] >> 12) << 14;
                    orank_row[pix - 2] |= (irank_row[cix] >> 10) << 14;
                    orank_row[pix - 3] |= (irank_row[cix] >> 8) << 14;
                    orank_row[pix - 4] |= (irank_row[cix] >> 6) << 14;
                    orank_row[pix - 5] |= (irank_row[cix] >> 4) << 14;
                    orank_row[pix - 6] |= (irank_row[cix] >> 2) << 14;
                    orank_row[pix - 7] |= (irank_row[cix] >> 0) << 14;

//                    {
//                        __v8hu v8 = (__v8hu)_mm_loadu_si128((__m128i *)&orank_row[pix - 7]);
//
//                        std::uint16_t hi = _mm_movemask_epi8((__m128i)v8);
//                        hi &= 0x2AAA;
//                        std::uint16_t lo = _mm_movemask_epi8(_mm_slli_epi16((__m128i)v8, 1));
//                        lo &= 0x2AAA;
//                        std::uint16_t unp = hi | (lo >> 1);
//                        assert(unp == irank_row[cix]);
//                    }
                }
                else
                {
                    orank_row.push_back(irank_row[cix]);
                }
            }

            for (auto cix = 0u; cix < NGENES; ++cix)
            {
                auto unp = [](std::uint16_t const * p, std::size_t ix) -> std::uint16_t
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
                auto real = irank_row[cix];
                auto ext = unp(orank_row.data(), cix);
                if (real != ext)
                {
                    std::cout << "real " << real << " ext " << ext << " pos " << cix << std::endl;
                }
                assert(real == ext);
            }

            fwrite(orank_row.data(), sizeof (orank_row[0]), orank_row.size(), orankf);
            orank_row.clear();
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    fclose(irankf);

    std::cout << "Done\n\n";

    std::cout << "Expected checksums:\n";
    std::cout << "sha1sum ranksBySigInvPacked\n\n";

    std::cout << "???  ranksBySigInvPacked\n\n";
}
