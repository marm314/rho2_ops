#!/bin/bash

name=`git rev-parse HEAD`
#echo $name

# cpp
echo '#include"gitver.h"' > gitver_new.cpp
echo " " >> gitver_new.cpp
echo "void gitversion(string &sha)" >> gitver_new.cpp
echo "{" >> gitver_new.cpp
echo " sha=\""$name"\";" >> gitver_new.cpp
echo "}" >> gitver_new.cpp

# Check if this gitver.cpp file is different (new SHA).
DIFF=$(diff gitver_new.cpp gitver.cpp)
if [ "$DIFF" != "" ] 
then
    echo "New SHA key"
    mv gitver_new.cpp gitver.cpp
fi
if [ "$DIFF" == "" ] 
then
    /bin/rm gitver_new.cpp
fi

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
