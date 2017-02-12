#define NDEBUG
#define NO_CMAPLIB_DEF
#define CREATE_INDEX_FILES

#pragma GCC target ("sse2")
#pragma GCC target ("sse4.2")

typedef unsigned int __v4su __attribute__ ((__vector_size__ (16)));

#include "query_parser.cpp"

#include "CMap2.hpp"
