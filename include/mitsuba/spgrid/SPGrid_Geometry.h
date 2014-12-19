//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class SPGrid_Geometry
//#####################################################################
#ifndef __SPGrid_Geometry_h__
#define __SPGrid_Geometry_h__

#include "std_array.h"

namespace SPGrid{

template<int dim> struct SPGrid_Geometry;

template<>
struct SPGrid_Geometry<3>
{
public:
    const unsigned xsize,ysize,zsize;                      // Dimensions requested
    const unsigned zsize_padded,ysize_padded,xsize_padded; // Dimensions allocated; adjusted for alignment and uniformity
    const unsigned block_xsize,block_ysize,block_zsize;    // Dimensions of data layout within block
    
    SPGrid_Geometry(const unsigned xsize_input,const unsigned ysize_input,const unsigned zsize_input,
        const unsigned block_xbits,const unsigned block_ybits,const unsigned block_zbits)
        :xsize(xsize_input),ysize(ysize_input),zsize(zsize_input),
        zsize_padded(Next_Power_Of_Two(max(xsize,ysize,zsize,1u<<block_zbits))),
        ysize_padded(max(Next_Power_Of_Two(max(ysize,xsize)),zsize_padded>>1,1u<<block_ybits)),
        xsize_padded(max(Next_Power_Of_Two(xsize),zsize_padded>>1,1u<<block_xbits)),
        block_xsize(1<<block_xbits),block_ysize(1<<block_ybits),block_zsize(1<<block_zbits)
    {}

    SPGrid_Geometry(SPGrid_Geometry &&other):
        xsize(other.xsize),
        ysize(other.ysize),
        zsize(other.zsize),
        zsize_padded(other.zsize_padded),
        ysize_padded(other.ysize_padded),
        xsize_padded(other.xsize_padded),
        block_xsize(other.block_xsize),
        block_ysize(other.block_ysize),
        block_zsize(other.block_zsize)
    {}
    
    unsigned long Padded_Volume() const
    {return (unsigned long)xsize_padded*(unsigned long)ysize_padded*(unsigned long)zsize_padded;}

    unsigned int Elements_Per_Block() const
    {return block_xsize*block_ysize*block_zsize;}

    std_array<unsigned,3> Block_Size() const
    {return std_array<unsigned,3>(block_xsize,block_ysize,block_zsize);}

    std_array<unsigned,3> Size() const
    {return std_array<unsigned,3>(xsize,ysize,zsize);}

//#####################################################################
    void Check_Bounds(const unsigned int i,const unsigned int j,const unsigned int k) const;
//#####################################################################
};

template<>
struct SPGrid_Geometry<2>
{
public:
    const unsigned xsize,ysize;               // Dimensions requested
    const unsigned ysize_padded,xsize_padded; // Dimensions allocated; adjusted for alignment and uniformity
    const unsigned block_xsize,block_ysize;   // Dimensions of data layout within block

    SPGrid_Geometry(const unsigned xsize_input,const unsigned ysize_input,const unsigned block_xbits,const unsigned block_ybits)
        :xsize(xsize_input),ysize(ysize_input),
        ysize_padded(Next_Power_Of_Two(max(xsize,ysize,1u<<block_ybits))),
        xsize_padded(max(Next_Power_Of_Two(xsize),ysize_padded>>1,1u<<block_xbits)),
        block_xsize(1<<block_xbits),block_ysize(1<<block_ybits)
    {}

    SPGrid_Geometry(SPGrid_Geometry &&other):
        xsize(other.xsize),
        ysize(other.ysize),
        ysize_padded(other.ysize_padded),
        xsize_padded(other.xsize_padded),
        block_xsize(other.block_xsize),
        block_ysize(other.block_ysize)
    {}
    
    unsigned long Padded_Volume() const
    {return (unsigned long)xsize_padded*(unsigned long)ysize_padded;}
    
    unsigned int Elements_Per_Block() const
    {return block_xsize*block_ysize;}

    std_array<unsigned,2> Block_Size() const
    {return std_array<unsigned,2>(block_xsize,block_ysize);}

    std_array<unsigned,2> Size() const
    {return std_array<unsigned,2>(xsize,ysize);}

//#####################################################################
    void Check_Bounds(const unsigned int i,const unsigned int j) const;
//#####################################################################
};
}

#endif
