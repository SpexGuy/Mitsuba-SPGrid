//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class HIERARCHICAL_RASTERIZER
//#####################################################################
#ifndef __HIERARCHICAL_RASTERIZER__
#define __HIERARCHICAL_RASTERIZER__

// #include <iostream>
// #include <utility>
// #include <stack>
// #include <list>
// #include <math.h>
// #include "std_array.h"

// using namespace std;

// namespace SPGrid{
// typedef std_array<int,3> Vec3;


template<class FUNCTOR>
class HIERARCHICAL_RASTERIZER
{
    FUNCTOR& functor;

public:

    HIERARCHICAL_RASTERIZER(FUNCTOR& functor_input)
        : functor(functor_input)
    {
    }

    void Iterate(GEOMETRY_BLOCK block)
    {
        int imin=block.imin.data[0];
        int jmin=block.imin.data[1];
        int kmin=block.imin.data[2];
        int imax=block.imax.data[0];
        int jmax=block.imax.data[1];
        int kmax=block.imax.data[2];

        float Xmin=block.Xmin.data[0];
        float Ymin=block.Xmin.data[1];
        float Zmin=block.Xmin.data[2];
        float Xmax=block.Xmax.data[0];
        float Ymax=block.Xmax.data[1];
        float Zmax=block.Xmax.data[2];

        if(!functor.Consume(block)) return;

        int half_size=(imax-imin)/2;
        float half_dx=(Xmax-Xmin)/2;

        std_array<int,3> root_i;
        std_array<float,3> root_X;

        root_i=std_array<int,3>  (imin          ,jmin          ,kmin          );
        root_X=std_array<float,3>(Xmin          ,Ymin          ,Zmin          );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));


        root_i=std_array<int,3>  (imin+half_size,jmin          ,kmin          );
        root_X=std_array<float,3>(Xmin+half_dx  ,Ymin          ,Zmin          );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin          ,jmin+half_size,kmin          );
        root_X=std_array<float,3>(Xmin          ,Ymin+half_dx  ,Zmin          );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin+half_size,jmin+half_size,kmin          );
        root_X=std_array<float,3>(Xmin+half_dx  ,Ymin+half_dx  ,Zmin          );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin          ,jmin          ,kmin+half_size);
        root_X=std_array<float,3>(Xmin          ,Ymin          ,Zmin+half_dx  );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin+half_size,jmin          ,kmin+half_size);
        root_X=std_array<float,3>(Xmin+half_dx  ,Ymin          ,Zmin+half_dx  );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin          ,jmin+half_size,kmin+half_size);
        root_X=std_array<float,3>(Xmin          ,Ymin+half_dx  ,Zmin+half_dx  );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));

        root_i=std_array<int,3>  (imin+half_size,jmin+half_size,kmin+half_size);
        root_X=std_array<float,3>(Xmin+half_dx  ,Ymin+half_dx  ,Zmin+half_dx  );
        Iterate(GEOMETRY_BLOCK(root_i,root_i+half_size,root_X,root_X+half_dx));


    }

//#####################################################################
};
#endif
