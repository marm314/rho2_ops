#!/bin/bash

name=`git rev-parse HEAD`
#echo $name

# cpp
echo '#include"gitver.h"' > gitver.cpp
echo " " >> gitver.cpp
echo "void gitversion(string &sha)" >> gitver.cpp
echo "{" >> gitver.cpp
echo " sha=\""$name"\";" >> gitver.cpp
echo "}" >> gitver.cpp

# header
echo '#ifndef _GIT_VERSION_H_' > gitver.h
echo '#define _GIT_VERSION_H_' >> gitver.h
echo " " >> gitver.h
echo '#include<iostream>' >> gitver.h
echo '#include"String_ops.h"' >> gitver.h
echo " " >> gitver.h
echo "using namespace std;" >> gitver.h
echo " " >> gitver.h
echo "void gitversion(string &sha);" >> gitver.h
echo "#endif // _GIT_VERSION_H_" >> gitver.h
