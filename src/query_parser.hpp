/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: query_parser.hpp
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
 * 2017-02-01   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_QUERY_PARSER_HPP_
#define SRC_QUERY_PARSER_HPP_

#include <cstdint>
#include <vector>
#include <string>

std::vector<std::uint32_t> parse_query(std::string const & s);

#endif /* SRC_QUERY_PARSER_HPP_ */
