/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: data_builder.cpp
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
 * 2017-02-14   wm              Initial version
 *
 ******************************************************************************/

#include <iostream>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <functional>

#undef NDEBUG

enum { NGENES = 10174 };
enum { NSIGS = 476251 };

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout << "I need two file names: scoresBySig and ranksBySig" << std::endl;
        return 1;
    }

    auto isigf = fopen(argv[1], "rb");
    if (isigf == nullptr)
    {
        std::cout << "Could not open for reading " << argv[1] << std::endl;
        return 1;
    }

    auto irankf = fopen(argv[2], "rb");
    if (irankf == nullptr)
    {
        std::cout << "Could not open for reading " << argv[2] << std::endl;
        return 1;
    }

    auto osigf = fopen("scoresBySigSortedAbsF32", "wb");
    if (osigf == nullptr)
    {
        std::cout << "Could not open for writing scoresBySigSortedAbsF32" << std::endl;
        return 1;
    }

    auto orankf = fopen("ranksBySigInv", "wb");
    if (orankf == nullptr)
    {
        std::cout << "Could not open for writing ranksBySigInv" << std::endl;
        return 1;
    }

    std::cout << "Assuming " << NSIGS << " signatures and " << NGENES << " genes" << std::endl;


    ////////////////////////////////////////////////////////////////////////////

    {
        double isig_row[NGENES];
        float osig_row[NGENES];

        for (auto six = 0u; six < NSIGS; ++six)
        {
            auto nread = fread(isig_row, sizeof (isig_row[0]), NGENES, isigf);
            assert(nread == NGENES);

            std::sort(std::begin(isig_row), std::end(isig_row), std::greater<double>());

            std::transform(std::begin(isig_row), std::end(isig_row), std::begin(osig_row),
                [](double x) -> float
                {
                    return std::abs(x);
                });
            fwrite(osig_row, sizeof (osig_row[0]), NGENES, osigf);
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    {
        int irank_row[NGENES];
        std::uint16_t orank_row[NGENES];

        for (auto six = 0u; six < NSIGS; ++six)
        {
            auto nread = fread(irank_row, sizeof (irank_row[0]), NGENES, irankf);
            assert(nread == NGENES);

            for (auto cix = 0u; cix < NGENES; ++cix)
            {
                irank_row[cix] = (irank_row[cix] << 16) | cix;
            }

            std::sort(std::begin(irank_row), std::end(irank_row));

            std::transform(std::begin(irank_row), std::end(irank_row), std::begin(orank_row),
                [](int x) -> std::uint16_t
                {
                    return x;
                });

            fwrite(orank_row, sizeof (orank_row[0]), NGENES, orankf);
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    fclose(isigf);
    fclose(irankf);

    std::cout << "Done\n\n";

    std::cout << "Expected checksums:\n";
    std::cout << "sha1sum scoresBySigSortedAbsF32\n\n";

    std::cout << "8b043a2f4af2be3e4c3d587565ed4b957c1a2dd5  scoresBySigSortedAbsF32\n\n";

    std::cout << "Expected checksums:\n";
    std::cout << "sha1sum ranksBySigInv\n\n";

    std::cout << "49c2d91d509e49ed86714561b4a93fee94cc39b7  ranksBySigInv\n\n";
}
