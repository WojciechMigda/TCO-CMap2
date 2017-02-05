/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: ground_truth.hpp
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

#ifndef SRC_GROUND_TRUTH_HPP_
#define SRC_GROUND_TRUTH_HPP_

#include <vector>
#include <string>

std::vector<double> load_gt_from_csv(std::vector<std::string> && vs);

#endif /* SRC_GROUND_TRUTH_HPP_ */
