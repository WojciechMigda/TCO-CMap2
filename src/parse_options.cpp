/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: parse_options.cpp
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

#include <boost/program_options.hpp>
#include <boost/any.hpp>

#include <iostream>

boost::program_options::variables_map
parse_options(int argc, char **argv)
{
    namespace po = boost::program_options;

    po::variables_map args;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("seed", po::value<int>()->default_value(1), "SEED for RNG")

        ("nsig", po::value<unsigned int>()->default_value(10), "number of signatures to match against")
        ("nskip", po::value<unsigned int>()->default_value(0), "number of signatures to skip")
//        ("up-query", po::value<std::string>()->default_value("../data/offline_query_up_n250.csv"), "Up query")
//        ("dn-query", po::value<std::string>()->default_value("../data/offline_query_down_n250.csv"), "Down query")
        //("gt-csv", po::value<std::string>()->default_value("../data/code/data/wtks_result_n10x10.csv"), "Ground truth")
        ("gt", po::value<std::string>()->default_value(""), "Ground truth")
        ("up-query", po::value<std::string>()->default_value("../data/code/data/query_up_n10.csv"), "Up query")
        ("dn-query", po::value<std::string>()->default_value("../data/code/data/query_down_n10.csv"), "Down query")
        ;
    try
    {
        po::store(po::parse_command_line(argc, argv, desc), args);
    }
    catch (po::error & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
    catch (boost::bad_any_cast & ex)
    {
        std::cerr << ex.what() << std::endl;
    }
    po::notify(args);

    if (args.count("help"))
    {
        std::cerr << desc << "\n";
        std::exit(1);
    }

    return args;
}
