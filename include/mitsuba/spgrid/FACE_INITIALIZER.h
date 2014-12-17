//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
// Class FACE_INITIALIZER
//#####################################################################
#ifndef __FACE_INITIALIZER__
#define __FACE_INITIALIZER__

#include "SPGrid_Set.h"
#include "SPGrid_Block_Iterator.h"
#include "FLUIDS_SIMULATION_FLAGS.h"

namespace SPGrid{

template<class T_STRUCT, int d>
class FACE_INITIALIZER
{
    typedef typename SPGrid_Allocator<T_STRUCT,d>::template Array<unsigned>::type Flag_array_type;
    typedef typename Flag_array_type::MASK Flag_array_mask;
    typedef SPGrid_Set<Flag_array_type> Set_type;

public:
    static void Flag_Active_Faces(Set_type& set)
    {
        int count = 0;
        static const int number_of_face_neighbors = (d==2) ? 4 : 6;
        unsigned long face_neighbor_shifts[number_of_face_neighbors];

        for (int axis = 0, face_neighbor = 0; axis < d; axis++) {
            for (int axis_shift = -1; axis_shift <= 1; axis_shift += 2) {
                std_array<int,d> shift;
                shift.data[axis]=axis_shift;
                face_neighbor_shifts[face_neighbor++]=Flag_array_mask::Linear_Offset(shift);
            }
        }
        set.Refresh_Block_Offsets();
        Flag_array_type flag_array = set.array;

        // Tag faces
        for(SPGrid_Block_Iterator<Flag_array_mask> iterator(set.Get_Blocks());iterator.Valid();iterator.Next())
            if(iterator.Data(flag_array) & SPGrid_Cell_Type_Interior) {
                count++;
                for(int face_neighbor=0;face_neighbor<number_of_face_neighbors;face_neighbor++)
                {
                    unsigned long neighbor_offset=Flag_array_mask::Packed_Add(iterator.Offset(),face_neighbor_shifts[face_neighbor]);
                    if(set.Is_Set(neighbor_offset,SPGrid_Cell_Type_Interior))
                        iterator.Data(flag_array) |= (SPGrid_Face_Minus_X_Active<<face_neighbor);
                }
            }

        //std::cout << "Processed " << count << " active cells.\n";
    }


//#####################################################################
};
}
#endif
