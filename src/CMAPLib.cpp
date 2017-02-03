/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: CMAPLib.cpp
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
 * 2017-02-03   wm              Initial version
 *
 ******************************************************************************/

#include "CMAPLib.hpp"
#include "likely.h"

#include <cstddef>
#include <vector>
#include <cstdio>
#include <cstdint>

using int_type = int;
using real_type = double;

template<typename Tp>
std::vector<Tp> loadFromFile(char const * fname, std::size_t pos, int nelem)
{
    std::vector<Tp> ret(nelem);

    FILE * ifile = fopen(fname, "rb");

    if (ifile)
    {
        fseek(ifile, pos, SEEK_SET);
        int nread = fread(&ret[0], sizeof (Tp), nelem, ifile);
        fclose(ifile);
        if (UNLIKELY(nread != nelem))
        {
            ret.resize(nread);
        }
    }

    return ret;
}

std::vector<double> loadFromDoubleFile(char const * fname, std::size_t pos, int nelem)
{
    return loadFromFile<double>(fname, pos, nelem);
}

std::vector<int_type> loadFromIntFile(char const * fname, std::size_t pos, int nelem)
{
    return loadFromFile<int_type>(fname, pos, nelem);
}

std::vector<double> CMAPLib::loadFromDoubleFile(std::string const & fname, std::size_t pos, int nelem)
{
    return ::loadFromDoubleFile(fname.c_str(), pos, nelem);
}

std::vector<int_type> CMAPLib::loadFromIntFile(std::string const & fname, std::size_t pos, int nelem)
{
    return ::loadFromIntFile(fname.c_str(), pos, nelem);
}


template<typename Tp>
int saveFile(char const * fname, std::vector<Tp> const & data)
{
    FILE * ofile = fopen(fname, "wb");

    if (ofile)
    {
        fwrite(data.data(), sizeof (Tp), data.size(), ofile);
        fclose(ofile);
    }

    return 0;
}

int saveIntFile(char const * fname, std::vector<int_type> const & data)
{
    return saveFile<int_type>(fname, data);
}

int saveDoubleFile(char const * fname, std::vector<double> const & data)
{
    return saveFile<double>(fname, data);
}

int CMAPLib::saveIntFile(std::string const & fname, std::vector<int_type> const & data)
{
    return ::saveIntFile(fname.c_str(), data);
}

int CMAPLib::saveDoubleFile(std::string const & fname, std::vector<double> const & data)
{
    return ::saveDoubleFile(fname.c_str(), data);
}
