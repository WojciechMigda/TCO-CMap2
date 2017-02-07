/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: ground_truth.cpp
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
 * 2017-02-05   wm              Initial version
 *
 ******************************************************************************/

#include "boost/spirit/include/qi_action.hpp"
#include "boost/spirit/include/qi_operator.hpp"

#include "boost/spirit/include/qi_parse.hpp"
#include "boost/spirit/include/qi_real.hpp"
#include "boost/spirit/include/qi_eps.hpp"
#include "boost/spirit/include/qi_char.hpp"

#include "cpplinq.hpp"

#include <cstddef>
#include <cmath>
#include <vector>

std::vector<double> load_gt_from_csv(std::vector<std::string> && vs, std::size_t nsig, std::size_t nskip)
{
    std::vector<double> ret;

    namespace q = ::cpplinq;
    q::from(vs)
        >> q::skip(1 + nskip)
        >> q::take(nsig)
        >> q::select([](std::string const & s)
            {
                return s.substr(s.find(',') + 1);
            })
        >> q::aggregate(&ret, [](std::vector<double> * seed, std::string const & s)
            {
                using namespace boost::spirit;
                auto first = s.c_str();
                auto last = first + s.size();

                bool ok = qi::phrase_parse(
                    first, last,
                    qi::double_ % ',',
                    qi::ascii::space,
                    *seed
                );
                return seed;
            })
        ;

    return ret;
}