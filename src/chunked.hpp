/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: chunked.hpp
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
 * 2017-02-10   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_CHUNKED_HPP_
#define SRC_CHUNKED_HPP_

#include "CMAPLib.hpp"

#include <cstddef>
#include <string>
#include <cstdint>
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>

std::string chunk_fname(std::string const & base, std::size_t num)
{
    return base + ".s" + std::to_string(num + 1);
}

void create_index_files(
    std::string const & ifile, std::string const & ofile_base,
    std::size_t NROWS, std::size_t NCOLS, std::size_t ROWS_PER_CHUNK)
{
    CMAPLib io;

    for (auto row = 0u; row < NROWS; row += ROWS_PER_CHUNK)
    {
        auto actual_rows = std::min(NROWS, row + ROWS_PER_CHUNK) - row;
        auto ichunk = io.loadFromIntFile(ifile, row * NCOLS, actual_rows * NCOLS);
        auto ibase = reinterpret_cast<std::uint32_t *>(ichunk.data());

        std::vector<int> ochunk(ichunk.size() / 2);
        auto obase = reinterpret_cast<std::uint16_t *>(ochunk.data());

        for (auto rix = 0u; rix < actual_rows; ++rix)
        {
            auto ibegin = ibase + (rix * NCOLS);
            auto obegin = obase + (rix * NCOLS);

            for (auto cix = 0u; cix < NCOLS; ++cix)
            {
                ibegin[cix] = (ibegin[cix] << 16) | cix;
            }

            std::sort(ibegin, ibegin + NCOLS);

            for (auto cix = 0u; cix < NCOLS; ++cix)
            {
                obegin[cix] = ibegin[cix];
            }
        }

        io.saveIntFile(chunk_fname(ofile_base, row / ROWS_PER_CHUNK), ochunk);
        std::cout << "Saved " << chunk_fname(ofile_base, row / ROWS_PER_CHUNK) << std::endl;
    }
}

void create_score_files(
    std::string const & ifile, std::string const & ofile_base,
    std::size_t NROWS, std::size_t NCOLS, std::size_t ROWS_PER_CHUNK)
{
    CMAPLib io;

    for (auto row = 0u; row < NROWS; row += ROWS_PER_CHUNK)
    {
        auto actual_rows = std::min(NROWS, row + ROWS_PER_CHUNK) - row;
        auto ichunk = io.loadFromDoubleFile(ifile, row * NCOLS, actual_rows * NCOLS);
        auto ibase = reinterpret_cast<double *>(ichunk.data());

        std::vector<int> ochunk(ichunk.size());
        auto obase = reinterpret_cast<float *>(ochunk.data());

        for (auto rix = 0u; rix < actual_rows; ++rix)
        {
            auto ibegin = ibase + (rix * NCOLS);
            auto obegin = obase + (rix * NCOLS);

            std::sort(ibegin, ibegin + NCOLS, std::greater<double>());

            for (auto cix = 0u; cix < NCOLS; ++cix)
            {
                obegin[cix] = std::abs(ibegin[cix]);
            }
        }

        io.saveIntFile(chunk_fname(ofile_base, row / ROWS_PER_CHUNK), ochunk);
        std::cout << "Saved " << chunk_fname(ofile_base, row / ROWS_PER_CHUNK) << std::endl;
    }
}

#endif /* SRC_CHUNKED_HPP_ */
