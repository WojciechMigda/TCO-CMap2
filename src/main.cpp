/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: main.cpp
 *
 * Description:
 *      TCO-CMap2
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

#include "CMap2.hpp"

#include "ground_truth.hpp"
#include "parse_options.hpp"

#include "cpplinq.hpp"

#include <vector>
#include <fstream>
#include <string>
#include <iostream>

std::vector<std::string>
read_file(const char * fname)
{
    std::ifstream fcsv(fname);
    std::vector<std::string> vcsv;

    for (std::string line; std::getline(fcsv, line); /* nop */)
    {
        vcsv.push_back(line);
    }
    fcsv.close();

    return vcsv;
}

int main(int argc, char **argv)
{
    auto const args = parse_options(argc, argv);

    auto query_selector = [](std::string const & s)
        {
            auto sec_comma = s.find(',', s.find(',') + 1);
            return s.substr(sec_comma + 1);
        };
    namespace q = ::cpplinq;
    auto up_query = q::from(read_file(args.at("up-query").as<std::string>().c_str()))
        >> q::select(query_selector)
        >> q::to_vector()
        ;

    auto dn_query = q::from(read_file(args.at("dn-query").as<std::string>().c_str()))
        >> q::select(query_selector)
        >> q::to_vector()
        ;

    ////////////////////////////////////////////////////////////////////////////


    std::cout << "[main] Loading ground-truth\n";
    auto const gt = args.at("gt").as<std::string>().size() > 0 ?
        load_gt_from_csv(
            read_file(args.at("gt").as<std::string>().c_str()),
            args.at("nsig").as<unsigned int>(),
            args.at("nskip").as<unsigned int>()) :
        std::vector<double>();
    std::cout << "[main] Done\n";

    CMAP2Updated solver(
        args.at("nsig").as<unsigned int>(),
        args.at("nskip").as<unsigned int>(),
        "../data/",
        gt.size() ? &gt : nullptr);

    auto genes = solver.m_cmap_lib.loadFromIntFile("../data/genes.bin", 0, 10174);
    solver.init(genes);

    solver.getWTKScomb(up_query, dn_query);

}
