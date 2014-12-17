//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __GEOMETRY_BLOCK__
#define __GEOMETRY_BLOCK__

#include "std_array.h"
using namespace SPGrid;

struct GEOMETRY_BLOCK
{
    std_array<int,3> imin;
    std_array<int,3> imax;
    std_array<float,3> Xmin;
    std_array<float,3> Xmax;

    GEOMETRY_BLOCK(std_array<int,3> imin_input,std_array<int,3> imax_input,std_array<float,3> Xmin_input,std_array<float,3> Xmax_input)
        :imin(imin_input),imax(imax_input),Xmin(Xmin_input),Xmax(Xmax_input) {}
};
#endif
