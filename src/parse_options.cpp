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
//        ("up-query", po::value<std::string>()->default_value("../data/offline_query_up_n250.csv"), "Up query")
//        ("dn-query", po::value<std::string>()->default_value("../data/offline_query_down_n250.csv"), "Down query")
        //("gt-csv", po::value<std::string>()->default_value("../data/code/data/wtks_result_n10x10.csv"), "Ground truth")
        ("gt-csv", po::value<std::string>()->default_value(""), "Ground truth")
        ("up-query", po::value<std::string>()->default_value("../data/code/data/query_up_n10.csv"), "Up query")
        ("dn-query", po::value<std::string>()->default_value("../data/code/data/query_down_n10.csv"), "Down query")


//        ("nfolds", po::value<int>()->default_value(3), "number of folds for cross-validation")
//
//        ("xgboost_params", po::value<bool>()->default_value(false), "XGBoost: override model params")
//        ("n_estimators", po::value<int>()->default_value(50), "XGBoost: number of estimators")
//        ("booster", po::value<std::string>()->default_value("gbtree"), "XGBoost: booster")
//        ("objective", po::value<std::string>()->default_value("rank:mmrf"), "XGBoost: objective")
//
//        ("colsample_bytree", po::value<float>()->default_value(1.0f),
//            "XGBoost: [gbtree] subsample ratio of columns when constructing each tree. range: (0,1].")
//        ("scale_pos_weight", po::value<float>()->default_value(1.0f),
//            "XGBoost: [gbtree] Control the balance of positive and negative weights, useful for unbalanced classes. A typical value to consider: sum(negative cases) / sum(positive cases).")
//        ("learning_rate", po::value<float>()->default_value(0.3f),
//            "XGBoost: [gbtree] step size shrinkage used in update to prevents overfitting. After each boosting step, we can directly get the weights of new features. and eta actually shrinks the feature weights to make the boosting process more conservative. range: [0,1]")
//        ("subsample", po::value<float>()->default_value(1.0f),
//            "XGBoost: [gbtree] subsample ratio of the training instance. Setting it to 0.5 means that XGBoost randomly collected half of the data instances to grow trees and this will prevent overfitting. range: (0,1].")
//        ("min_child_weight", po::value<float>()->default_value(1.0f),
//            "XGBoost: [gbtree] minimum sum of instance weight (hessian) needed in a child. If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. In linear regression mode, this simply corresponds to minimum number of instances needed to be in each node. The larger, the more conservative the algorithm will be. range: [0,∞]")
//        ("max_depth", po::value<int>()->default_value(6),
//            "XGBoost: [gbtree] maximum depth of a tree, increase this value will make the model more complex / likely to be overfitting. range: [1,∞]")
//        ("num_pairsample", po::value<int>()->default_value(1), "XGBoost: num_pairsample")
//
//        ("reg_lambda", po::value<float>()->default_value(0.0f),
//            "XGBoost: [gblinear] L2 regularization term on weights, increase this value will make model more conservative.")
//        ("reg_alpha", po::value<float>()->default_value(0.0f),
//            "XGBoost: [gblinear] L1 regularization term on weights, increase this value will make model more conservative.")
//        ("reg_lambda_bias", po::value<float>()->default_value(0.0f),
//            "XGBoost: [gblinear] L2 regularization term on bias")
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
