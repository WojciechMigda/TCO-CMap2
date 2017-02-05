#!/bin/sh

/repo/topcoder/Amalgamate/Amalgamate \
-w "*.cpp;*.c;*.cc;*.h;*.hpp" \
-i . \
-i ../cpplinq/CppLinq/ \
amalgamate.cpp submission.cpp

sed -i  "s/_ZStL8__ioinit/_ZStL8__ioinit_dummy/g" submission.cpp

#g++ -std=c++11 -c submission.cpp
g++ -std=c++11 -c wrapped_submission.cpp -W -Wall -Wno-sign-compare -O2 -s -pipe -mmmx -msse -msse2 -msse3 -march=native -o submission.o

#g++ -std=c++11 -c submission.cpp -W -Wall -Wno-sign-compare -Os -s -pipe #-mmmx -msse -msse2 -msse3
