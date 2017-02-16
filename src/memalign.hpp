/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: memalign.hpp
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
 * 2017-02-16   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_MEMALIGN_HPP_
#define SRC_MEMALIGN_HPP_

#include "likely.h"

#include <cstdlib>
#include <cstddef>

static inline
void * aligned_malloc(std::size_t a, std::size_t sz)
{
    void * p;
    auto ret = posix_memalign(&p, a, sz);
    if (LIKELY(ret == 0))
    {
        return p;
    }
    else
    {
        return nullptr;
    }
}


#endif /* SRC_MEMALIGN_HPP_ */
