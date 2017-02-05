/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: CMAPLib.hpp
 *
 * Description:
 *      Interface to score/rank file I/O
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

#ifndef SRC_CMAPLIB_HPP_
#define SRC_CMAPLIB_HPP_

#include <vector>
#include <cstddef>
#include <string>

#ifndef NO_CMAPLIB_DEF
class CMAPLib
{
public:
    std::vector<double>
    loadFromDoubleFile(std::string const & fname, std::size_t pos, int nelem);

    std::vector<int>
    loadFromIntFile(std::string const & fname, std::size_t pos, int nelem);

    int
    saveIntFile(std::string const & fname, std::vector<int> const & data);

    int
    saveDoubleFile(std::string const & fname, std::vector<double> const & data);
};
#endif

#endif /* SRC_CMAPLIB_HPP_ */
