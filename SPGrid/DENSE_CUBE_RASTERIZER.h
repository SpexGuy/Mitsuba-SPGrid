//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class DENSE_CUBE_RASTERIZER
//#####################################################################
#ifndef __DENSE_CUBE_RASTERIZER__
#define __DENSE_CUBE_RASTERIZER__
#include <stdio.h>
#include <iostream>

#include "std_array.h"

using namespace SPGrid;

template<class T_SET>
struct DENSE_CUBE_RASTERIZER
{
    static int Rasterize(T_SET& set, int size) {
        int count = 0;
        for (int x = 2; x < size - 2; ++x)
        for (int y = 2; y < size - 2; ++y)
        for (int z = 2; z < size - 2; ++z) {
            set.Mask(std_array<int,3>(x,y,z), SPGrid_Cell_Type_Interior);
            count++;
        }
        return count;
    }
//#####################################################################
};
#endif
