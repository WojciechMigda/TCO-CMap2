/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: CMap2.hpp
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
 * 2017-01-31   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_CMAP2_HPP_
#define SRC_CMAP2_HPP_

#include "CMAPLib.hpp"
#include "query_parser.hpp"

#include "cpplinq.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cstddef>

#undef NSIG

using size_type = std::size_t;

//#define TESTING
//
//#ifdef TESTING
//#include "array2d.hpp"
//#include <utility>
//#include <valarray>
//using real_type = double;
//using array_type = num::array2d<real_type>;
//using varray_type = std::valarray<real_type>;
//
//std::vector<std::string>
//read_file(const char * fname);
//#endif

#include <sys/time.h>
long int timestamp()
{
    timeval tv;
    gettimeofday(&tv, NULL);

    return tv.tv_sec;// + tv.tv_usec / 1e6;
}

struct CMAP2Updated
{
    enum {NGENES = 10174};

    CMAP2Updated(size_type nsig = 476251, std::string path = "", std::vector<double> const * gt = nullptr) :
        NSIG(nsig), PATH(path), m_gt(gt)
    {}

    int init(std::vector<int> const & genes);
    std::vector<double> getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);
    std::vector<double> getWTKScomb_(std::vector<std::string> & q_up, std::vector<std::string> & q_dn);

    size_type const NSIG;
    std::string const PATH;
    std::vector<double> const * m_gt;
    std::unordered_map<int, int> m_gene_to_idx;
    CMAPLib m_cmap_lib;
};

int
CMAP2Updated::init(std::vector<int> const & genes)
{
    m_gene_to_idx.clear();
    namespace q = ::cpplinq;
    q::from(genes)
        >> q::zip_with(q::detail::int_range(0, INT_MAX))
        >> q::aggregate(&m_gene_to_idx,
            [](std::unordered_map<int, int> * seed, std::pair<int, int> const & kv)
            {
                seed->emplace(kv.first, kv.second);
                return seed;
            })
        ;

    // cache warm up ??? TODO

    return 0;
}


std::vector<double>
CMAP2Updated::getWTKScomb(std::vector<std::string> & q_up, std::vector<std::string> & q_dn)
{
    const auto time0 = timestamp();

    assert(q_up.size() == q_dn.size());
    auto const NQRY = q_up.size();

    auto q_indexer = [this](std::vector<std::string> & qry)
        {
            namespace q = ::cpplinq;
            return q::from(qry)
                >> q::select([this](std::string const & s)
                    {
                        auto q_parsed = parse_query(s);
                        std::transform(q_parsed.begin(), q_parsed.end(), q_parsed.begin(),
                            [this](int gene_id){ return this->m_gene_to_idx.at(gene_id); });
                        return q_parsed;
                    })
                >> q::to_vector();
        };

    auto const q_up_indexed = q_indexer(q_up);
    auto const q_dn_indexed = q_indexer(q_dn);

    std::vector<double> sig;
    std::vector<int> ranks;

    std::vector<double> ret;
    ret.reserve(q_up.size() * NSIG);

    for (auto six = 0u; six < NSIG; ++six)
    {
        sig = m_cmap_lib.loadFromDoubleFile(PATH + "scoresBySig", six * NGENES, NGENES);
        ranks = m_cmap_lib.loadFromIntFile(PATH + "ranksBySig", six * NGENES, NGENES);
        // cache i/o ???

        for (auto qix = 0u; qix < NQRY; ++qix)
        {
            auto & q_up_indices = q_up_indexed[qix];
            auto & q_dn_indices = q_dn_indexed[qix];

            auto const QUSIZE = q_up_indices.size();
            auto const QDSIZE = q_dn_indices.size();

            std::vector<std::pair<int, double>> up_ranks;
            up_ranks.reserve(QUSIZE);
            double Sum_up_abs_scores = 0.;
            std::transform(q_up_indices.cbegin(), q_up_indices.cend(), std::back_inserter(up_ranks),
                [&ranks, &sig, &Sum_up_abs_scores](int ix)
                {
                    auto abs_score = std::abs(sig.at(ix));
                    Sum_up_abs_scores += abs_score;
                    return std::make_pair(ranks.at(ix), abs_score);
                });

            std::vector<std::pair<int, double>> dn_ranks;
            dn_ranks.reserve(QDSIZE);
            double Sum_dn_abs_scores = 0.;
            std::transform(q_dn_indices.cbegin(), q_dn_indices.cend(), std::back_inserter(dn_ranks),
                [&ranks, &sig, &Sum_dn_abs_scores](int ix)
                {
                    auto abs_score = std::abs(sig.at(ix));
                    Sum_dn_abs_scores += abs_score;
                    return std::make_pair(ranks.at(ix), abs_score);
                });

            std::qsort(&up_ranks[0], up_ranks.size(), sizeof (up_ranks.front()),
                [](const void * avp, const void * bvp)
                {
                    auto ap = static_cast<std::pair<int, double> const *>(avp);
                    auto bp = static_cast<std::pair<int, double> const *>(bvp);
                    return (ap->first > bp->first) - (ap->first < bp->first);
                }
            );
            std::qsort(&dn_ranks[0], dn_ranks.size(), sizeof (dn_ranks.front()),
                [](const void * avp, const void * bvp)
                {
                    auto ap = static_cast<std::pair<int, double> const *>(avp);
                    auto bp = static_cast<std::pair<int, double> const *>(bvp);
                    return (ap->first > bp->first) - (ap->first < bp->first);
                }
            );

            double const upenalty = -1. / (NGENES - QUSIZE);
            double wkts_up = 0.;
            {
                //volatile
                double _acc = 0.;
                double _extr = 0.;
                int prev = 0; // rank indices are 1-based
                for (auto const & rs : up_ranks)
                {
                    _acc += (rs.first - prev - 1) * upenalty;
                    //asm("popcnt %rax,%rax");
//                    for (auto aix = 0u; aix < (rs.first - prev - 1); ++aix)
//                    {
//                        _acc += penalty;
//                    }
                    auto abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_up = _acc;
                    }
                    prev = rs.first;

                    _acc += rs.second / Sum_up_abs_scores;
                    abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_up = _acc;
                    }
                }
//                if (prev != NGENES)
//                {
//                    _acc += (NGENES - prev) * penalty;
//                    auto const abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_up = _acc;
//                    }
//                }
            }

            double const dpenalty = -1. / (NGENES - QDSIZE);
            double wkts_dn = 0.;
            {
                double _acc = 0.;
                double _extr = 0.;
                int prev = 0; // rank indices are 1-based
                for (auto const & rs : dn_ranks)
                {
                    _acc += (rs.first - prev - 1) * dpenalty;
                    auto abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_dn = _acc;
                    }
                    prev = rs.first;

                    _acc += rs.second / Sum_dn_abs_scores;
                    abs_acc = std::abs(_acc);
                    if (abs_acc > _extr)
                    {
                        _extr = abs_acc;
                        wkts_dn = _acc;
                    }
                }
//                if (prev != NGENES)
//                {
//                    _acc += (NGENES - prev) * penalty;
//                    auto const abs_acc = std::abs(_acc);
//                    if (abs_acc > _extr)
//                    {
//                        _extr = abs_acc;
//                        wkts_dn = _acc;
//                    }
//                }
            }

            auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;
            //std::cout << "wkts : " << wkts << "\n";
            ret.push_back(wkts);
            if (m_gt)
            {
                auto ref = (*m_gt)[six * NQRY + qix];
                if (std::abs(ref - wkts) >= 0.001)
                {
                    std::cout << "Ref: " << ref << " WKTS: " << wkts << " up/dn " << wkts_up << ' ' << wkts_dn << '\n';
                }
            }
        }
    }

    return ret;
}

////////////////////////////////////////////////////////////////////////////////

//std::vector<double>
//CMAP2::getWTKScomb_(std::vector<std::string> & q_up, std::vector<std::string> & q_dn)
//{
//    constexpr auto NCOL = 10;
//    constexpr auto NROW = 10174;
//    using rank_type = std::size_t;
//
//#ifdef TESTING
//    auto scores_vs = read_file("../data/code/data/score_n10x10174.csv");
//    namespace q = ::cpplinq;
//    auto score_headers_str = scores_vs.front().substr(scores_vs.front().find(',') + 1);
//    std::cout << "score_headers_str : " << score_headers_str << "\n";
//    std::unordered_map<int, int> gene_to_idx;
//    q::from(scores_vs)
//        >> q::skip(1)
//        >> q::select([](std::string const & s){ return s.substr(0, s.find(',')); })
//        >> q::zip_with(q::detail::int_range(0, INT_MAX))
//        >> q::aggregate(&gene_to_idx, [](std::unordered_map<int, int> * seed, std::pair<std::string, int> const & kv){ seed->emplace(std::stoi(kv.first), kv.second); return seed; })
//        ;
//    std::cout << "Gene '100' index shall be 0 : " << gene_to_idx.at(100) << "\n";
//    auto string_as_varray = [NCOL](std::string const & s)
//        {
//            real_type row[NCOL];
//            sscanf(s.c_str(), "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", row + 0, row + 1, row + 2, row + 3, row + 4, row + 5, row + 6, row + 7, row + 8, row + 9);
//            return varray_type(row, NCOL);
//        };
//    auto varray_reducer = [](array_type * seed, std::pair<varray_type, int> const & zipped)
//        {
//            (*seed)[seed->row(zipped.second)] = zipped.first;
//            return seed;
//        };
//    array_type scores(num::shape_type(NROW, NCOL), 0.f);
//    q::from(scores_vs)
//        >> q::skip(1)
//        >> q::select([](std::string const & s){ return s.substr(s.find(',') + 1); })
//        >> q::select(string_as_varray)
//        >> q::zip_with(q::detail::int_range(0, INT_MAX))
//        >> q::aggregate(&scores, varray_reducer)
//        ;
//    std::cout << "first row of scores : ";
//    varray_type row0 = scores[scores.row(0)];
//    std::copy(std::begin(row0), std::end(row0), std::ostream_iterator<real_type>(std::cout, ", "));
//    std::cout << '\n';
//
//    num::array2d<rank_type> ranks(num::shape_type(NROW, NCOL));
//    for (size_type ix = 0; ix < NCOL; ++ix)
//    {
//        varray_type col = scores[scores.column(ix)];
//        auto r = q::from(std::vector<double>(std::begin(col), std::end(col)))
//            >> q::zip_with(q::detail::int_range(0, INT_MAX))
//            >> q::orderby_descending([](std::pair<double, int> const & vi){ return vi.first; })
//            >> q::zip_with(q::detail::int_range(1, INT_MAX))
//            >> q::orderby_ascending([](std::pair<std::pair<double, int>, int> const & vii){ return vii.first.second; })
//            >> q::select([](std::pair<std::pair<double, int>, int> const & vii) -> rank_type { return vii.second; })
//            >> q::to_vector(NROW)
//            ;
//        ranks[ranks.column(ix)] = std::valarray<rank_type>(r.data(), r.size());
//    }
//
//#endif
//
//    auto q_up_parsed = parse_query(q_up.front());
//    auto q_up_indices = q::from(q_up_parsed)
//        >> q::select([&gene_to_idx](std::uint32_t gene) -> size_type { return gene_to_idx.at(gene); })
//        >> q::to_vector(q_up_parsed.size())
//        ;
//
//    std::valarray<rank_type> up_ranks = static_cast<std::valarray<rank_type>>(ranks[ranks.column(0)])[std::valarray<size_type>(q_up_indices.data(), q_up_indices.size())];
//    std::valarray<double> run_stats(-1. / (NROW - q_up_parsed.size()), NROW);
//    std::valarray<double> up_abs_scores = std::abs<double>(static_cast<std::valarray<double>>(scores[scores.column(0)])[std::valarray<size_type>(q_up_indices.data(), q_up_indices.size())]);
//    auto const Sum_up_abs_scores = up_abs_scores.sum();
//    run_stats[up_ranks] = up_abs_scores / Sum_up_abs_scores;
//
//    double _acc = 0.;
//    double _extr = 0.;
//    double wkts_up = 0.;
//    auto it = std::begin(run_stats);
//    _acc = *it++;
//    wkts_up = _acc;
//    _extr = std::abs(_acc);
//    for (; it != std::end(run_stats); ++it)
//    {
//        _acc += *it;
//        auto const abs_acc = std::abs(_acc);
//        if (abs_acc > _extr)
//        {
//            _extr = abs_acc;
//            wkts_up = _acc;
//        }
//    }
//    std::cout << "wkts_up : " << wkts_up << "\n";
//
//
//    auto q_dn_parsed = parse_query(q_dn.front());
//    auto q_dn_indices = q::from(q_dn_parsed)
//        >> q::select([&gene_to_idx](std::uint32_t gene) -> size_type { return gene_to_idx.at(gene); })
//        >> q::to_vector(q_dn_parsed.size())
//        ;
//
//    std::valarray<rank_type> dn_ranks = static_cast<std::valarray<rank_type>>(ranks[ranks.column(0)])[std::valarray<size_type>(q_dn_indices.data(), q_dn_indices.size())];
//    run_stats = -1. / (NROW - q_dn_parsed.size());
//    std::valarray<double> dn_abs_scores = std::abs<double>(static_cast<std::valarray<double>>(scores[scores.column(0)])[std::valarray<size_type>(q_dn_indices.data(), q_dn_indices.size())]);
//    auto const Sum_dn_abs_scores = dn_abs_scores.sum();
//    run_stats[dn_ranks] = dn_abs_scores / Sum_dn_abs_scores;
//
//    _acc = 0.;
//    _extr = 0.;
//    double wkts_dn = 0.;
//    it = std::begin(run_stats);
//    _acc = *it++;
//    wkts_dn = _acc;
//    _extr = std::abs(_acc);
//    for (; it != std::end(run_stats); ++it)
//    {
//        _acc += *it;
//        auto const abs_acc = std::abs(_acc);
//        if (abs_acc > _extr)
//        {
//            _extr = abs_acc;
//            wkts_dn = _acc;
//        }
//    }
//    std::cout << "wkts_dn : " << wkts_dn << "\n";
//
//    auto wkts = (wkts_dn * wkts_up < 0.) ? (wkts_up - wkts_dn) / 2 : 0.;
//    std::cout << "wkts : " << wkts << "\n";
//
//    return std::vector<double>();
//}

#endif /* SRC_CMAP2_HPP_ */
