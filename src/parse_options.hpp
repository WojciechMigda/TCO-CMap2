/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: parse_options.hpp
 *
 * Description:
 *      Parse options using boost.program_options
 *
 * Authors:
 *          Wojciech Migda (wm)
 *
 *******************************************************************************
 * History:
 * --------
 * Date         Who  Ticket     Description
 * ----------   ---  ---------  ------------------------------------------------
 * 2017-01-30   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_PARSE_OPTIONS_HPP_
#define SRC_PARSE_OPTIONS_HPP_

#include <boost/program_options/variables_map.hpp>

boost::program_options::variables_map
parse_options(int argc, char **argv);


#endif /* SRC_PARSE_OPTIONS_HPP_ */
