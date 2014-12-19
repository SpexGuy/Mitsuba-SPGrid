#ifndef __Save_Helper__
#define __Save_Helper__

#include <iostream>
#include <fstream>
#include "SPGrid_Utilities.h"
#include "SPGrid_Allocator.h"

namespace SPGrid{

template<class T, unsigned int dim> // TODO: enforce SET_T comes from ALLOC_T, dim == alloc_t.dim
struct Save_Helper
{
    typedef SPGrid_Allocator<T, dim> ALLOC_T;
public:

    template<class SET_T>
    static bool Save_Mask(const ALLOC_T &alloc, SET_T &set, const std::string &filename) {
        std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
        if (!out.is_open())
            return false;
        // first line is the grid size
        std_array<unsigned,dim> size = alloc.Size();
        Static_Assert(dim > 0);
        out << size(0);
        for (int c = 1; c < dim; c++)
            out << ' ' << size(c);
        out << std::endl;
        // then the offsets
        set.Refresh_Block_Offsets();
        const unsigned long *blocks = set.Get_Blocks().first;
        unsigned numBlocks = set.Get_Blocks().second;
        for (int c = 0; c < numBlocks; c++) {
            out << blocks[c] << std::endl;
        }
        out.close();
        return true;
    }

    template<class T_FIELD, class T1, class SET_T> // TODO: somehow ensure that the set and array come from the same allocator type
    static bool Save_Data(const ALLOC_T &alloc, T_FIELD T1::* field, SET_T &set, const std::string &filename) {
        typedef  typename EnableIfSame<typename ALLOC_T::template Array<const T_FIELD>::type, T1, T>::type  Array_type;
        std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
        if (!out.is_open())
            return false;
        set.Refresh_Block_Offsets();
        const unsigned long *blocks = set.Get_Blocks().first;
        const unsigned numBlocks = set.Get_Blocks().second;
        const unsigned elementsPerBlock = alloc.Elements_Per_Block();
        Array_type data = alloc.Get_Array(field);
        const T_FIELD *rawData = reinterpret_cast<const T_FIELD*>(data.Get_Data_Ptr());
        for (int c = 0; c < numBlocks; c++) {
            out.write(reinterpret_cast<const char*>((unsigned long)rawData + blocks[c]), elementsPerBlock * sizeof(T_FIELD));
        }
        return true;
    }
};

}
#endif
