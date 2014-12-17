//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class ADAPTIVE_SPHERE_RASTERIZER
//#####################################################################
#ifndef __ADAPTIVE_SPHERE_RASTERIZER__
#define __ADAPTIVE_SPHERE_RASTERIZER__
#include <stdio.h>
#include <iostream>

#include "FLUIDS_SIMULATION_FLAGS.h"
#include "SPGrid_Set.h"
#include "GEOMETRY_BLOCK.h"
#include "std_array.h"

using namespace SPGrid;

template<class T_SET>
struct ADAPTIVE_SPHERE_RASTERIZER
{
    std_array<float,3> center;
    float inner_radius;
    float outer_radius;
    float inner_radius_squared;
    float outer_radius_squared;
    int total_active;
    T_SET& set;

public:
    ADAPTIVE_SPHERE_RASTERIZER(T_SET& set_input, std_array<float,3> center_input,float inner_radius_input,float outer_radius_input)
        :set(set_input), center(center_input),inner_radius(inner_radius_input),outer_radius(outer_radius_input),
         inner_radius_squared(inner_radius*inner_radius),outer_radius_squared(outer_radius*outer_radius), total_active(0) {}

    static float Magnitude_Squared(const std_array<float,3> v){
        return v.data[0]*v.data[0]+v.data[1]*v.data[1]+v.data[2]*v.data[2];
    }

    static float Clamp(const float x,const float xmin,const float xmax)
    {if(x<=xmin) return xmin;else if(x>=xmax) return xmax;else return x;}

    static std_array<float,3> Clamp(const std_array<float,3> v,const std_array<float,3> vmin,const std_array<float,3> vmax)
    {return std_array<float,3>(Clamp(v.data[0],vmin.data[0],vmax.data[0]),Clamp(v.data[1],vmin.data[1],vmax.data[1]),Clamp(v.data[2],vmin.data[2],vmax.data[2]));}

    float Maximum_Distance_Squared(const std_array<float,3> Xmin,const std_array<float,3> Xmax)
    {return Maximum_Distance_Squared(Xmin.data[0],Xmax.data[0],Xmin.data[1],Xmax.data[1],Xmin.data[2],Xmax.data[2]);}

    float Maximum_Distance_Squared(float Xmin,float Xmax,float Ymin,float Ymax,float Zmin,float Zmax)
    {
        float distance_squared=0.f;
        std_array<float,3> dX;
        float new_distance_squared;

        dX=std_array<float,3>(Xmin,Ymin,Zmin)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;        

        dX=std_array<float,3>(Xmax,Ymin,Zmin)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmin,Ymax,Zmin)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmax,Ymax,Zmin)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmin,Ymin,Zmax)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmax,Ymin,Zmax)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmin,Ymax,Zmax)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        dX=std_array<float,3>(Xmax,Ymax,Zmax)-center;
        new_distance_squared=Magnitude_Squared(dX);
        if(new_distance_squared>distance_squared) distance_squared=new_distance_squared;

        return distance_squared;
    }

    float Minimum_Distance_Squared(const std_array<float,3> Xmin,const std_array<float,3> Xmax)
    {
        return Magnitude_Squared(Clamp(center,Xmin,Xmax)-center);
    }

    bool Consume(GEOMETRY_BLOCK block)
    {
        std_array<int,3> lengths = block.imax-block.imin;
        if (lengths == std_array<int,3>(1)) {
            set.Mask(block.imax, SPGrid_Cell_Type_Interior);
            total_active++;
            return false;
        }

        float min_dist_sq = Minimum_Distance_Squared(block.Xmin, block.Xmax);
        if (min_dist_sq > outer_radius_squared) return false;  // outside band, don't recurse

        float max_dist_sq = Maximum_Distance_Squared(block.Xmin, block.Xmax);
        if (max_dist_sq < inner_radius_squared) return false;  // inside band, don't recurse

        return true;
        
    }
//#####################################################################
};
#endif
