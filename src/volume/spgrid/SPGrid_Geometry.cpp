//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class SPGrid_Geometry
//#####################################################################
#include <mitsuba/spgrid/SPGrid_Geometry.h>

using namespace SPGrid;
//#####################################################################
// Function Check_Bounds
//#####################################################################
void SPGrid_Geometry<3>::
Check_Bounds(const unsigned int i,const unsigned int j,const unsigned int k) const
{
    if(i>=xsize_padded || j>=ysize_padded || k>=zsize_padded)
        FATAL_ERROR("Array indices ("+Value_To_String(i)+","+Value_To_String(j)+","+Value_To_String(k)+") out of hard bounds");
    if(i>=xsize || j>=ysize || k>=zsize)
        FATAL_ERROR("Array indices ("+Value_To_String(i)+","+Value_To_String(j)+","+Value_To_String(k)+") out of soft bounds");
}
//#####################################################################
void SPGrid_Geometry<2>::
Check_Bounds(const unsigned int i,const unsigned int j) const
{
    if(i>=xsize_padded || j>=ysize_padded)
        FATAL_ERROR("Array indices ("+Value_To_String(i)+","+Value_To_String(j)+") out of hard bounds");
    if(i>=xsize || j>=ysize)
        FATAL_ERROR("Array indices ("+Value_To_String(i)+","+Value_To_String(j)+") out of soft bounds");
}
//#####################################################################
