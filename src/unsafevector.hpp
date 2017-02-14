/*******************************************************************************
 * Copyright (c) 2017 Wojciech Migda
 * All rights reserved
 * Distributed under the terms of the MIT License
 *******************************************************************************
 *
 * Filename: qvector.hpp
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
 * 2017-02-14   wm              Initial version
 *
 ******************************************************************************/

#ifndef SRC_UNSAFEVECTOR_HPP_
#define SRC_UNSAFEVECTOR_HPP_

#include "likely.h"

#include <cstdint>
#include <cstddef>
#include <cstdlib>

// unsafe, because it doesn't check bounds
template<typename Tp>
struct unsafe_vector
{
    using value_type = Tp;
    using pointer_type = value_type *;
    using size_type = std::size_t;
    using const_reference = value_type const &;

    unsafe_vector()
    {
        m_sz = 0;
        m_cap = 200;
        m_v = (pointer_type)malloc(m_cap * sizeof (value_type));
    }

    ~unsafe_vector()
    {
        free(m_v);
    }

    const_reference operator[](size_type pos) const
    {
        return m_v[pos];
    }

    void push_back(const std::uint16_t & value)
    {
        if (UNLIKELY(m_sz == m_cap))
        {
            m_cap *= 2;
            m_v = (pointer_type)realloc(m_v, m_cap * sizeof (value_type));
        }
        m_v[m_sz++] = value;
    }

    void push_back(value_type && value)
    {
        if (UNLIKELY(m_sz == m_cap))
        {
            m_cap *= 2;
            m_v = (pointer_type)realloc(m_v, m_cap * sizeof (value_type));
        }
        m_v[m_sz++] = value;
    }

    void reserve(size_type new_cap)
    {
        if (new_cap > m_cap)
        {
            m_cap = new_cap;
            m_v = (pointer_type)realloc(m_v, m_cap * sizeof (value_type));
        }
    }

    size_type size() const
    {
        return m_sz;
    }

    void clear()
    {
        m_sz = 0;
    }

    std::uint32_t m_sz;
    std::uint32_t m_cap;
    pointer_type m_v;
};

#endif /* SRC_UNSAFEVECTOR_HPP_ */
