/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: gt_check.cpp
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
 * 2017-02-16   wm              Initial version
 *
 ******************************************************************************/

#include <iostream>
#include <cstdio>
#include <cmath>

enum { NSIGS = 476251 };

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout << "I need two ground truth FP32 binary files" << std::endl;
        return 1;
    }

    std::cout << "Assuming " << NSIGS << " signatures" << std::endl;

    FILE * if1 = fopen(argv[1], "rb");
    if (if1 == nullptr)
    {
        std::cout << "Failed to open " << argv[1] << " for reading" << std::endl;
        return 1;
    }

    FILE * if2 = fopen(argv[2], "rb");
    if (if2 == nullptr)
    {
        std::cout << "Failed to open " << argv[2] << " for reading" << std::endl;
        return 1;
    }

    fseek(if1, 0, SEEK_END);
    fseek(if2, 0, SEEK_END);

    auto const if1_sz = ftell(if1);
    auto const if2_sz = ftell(if2);

    if (if1_sz != if2_sz)
    {
        std::cout << "Size mismatch between files " << argv[1] << " and " << argv[2] << std::endl;
        return 1;
    }
    rewind(if1);
    rewind(if2);

    if ((if1_sz % sizeof (float) != 0 ) || ((if1_sz / sizeof (float)) % NSIGS) != 0)
    {
        std::cout << "Number of FP32 records in file should be a multiply of number of signatures" << std::endl;
        return 1;
    }

    auto const NQRY = if1_sz / NSIGS;
    float b1[NSIGS];
    float b2[NSIGS];
    for (auto six = 0u; six < if1_sz / NQRY; ++six)
    {
        auto const nread1 = fread(b1, sizeof (b1[0]), NQRY, if1);
        auto const nread2 = fread(b2, sizeof (b2[0]), NQRY, if2);

        if (nread1 != nread2)
        {
            std::cout << "Read error" << std::endl;
        }

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            if (std::abs(b1[qix] - b2[qix]) >= 1e-3)
            {
                std::cout << "Sig: " << six + 1 << ", Query: " << qix + 1 \
                    << " Ref 1: " << b1[qix] << " Ref 2: " << b2[qix] << std::endl;
            }
        }
    }

    fclose(if1);
    fclose(if2);
}
