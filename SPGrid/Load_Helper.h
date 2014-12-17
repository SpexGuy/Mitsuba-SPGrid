#ifndef __Load_Helper__
#define __Load_Helper__

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "SPGrid_Allocator.h"
#include "SPGrid_Set.h"

using std::ios;
using std::ifstream;

namespace SPGrid{

template<class T, unsigned dim>
struct Load_Helper
{
    typedef SPGrid_Allocator<T, dim> Allocator_type;
public:
    Allocator_type &alloc;

    Load_Helper(Allocator_type &alloc_in)
        : alloc(alloc_in)
    {}

    template<class SET_T>
    bool Load_Mask(SET_T &set, const std::string &filename) {
        ifstream in(filename.c_str(), ios::in);
        if (!in.is_open())
            return false;

        // first line is the block size - verify that it is the same
        std::string dims;
        std_array<unsigned, dim> size = alloc.Block_Size();
        if (!std::getline(in, dims)) return false;
        std::istringstream firstLine(dims);
        unsigned length;
        unsigned dimCount = 0;
        while(firstLine >> length) {
            if (dimCount >= dim)
                return false;
            if (length != size(dimCount))
                return false;
            ++dimCount;
        }
        if (dimCount != dim)
            return false;

        // then the offsets
        unsigned long offset;
        std::vector<unsigned long> block_offsets;
        while (in >> offset) {
            block_offsets.push_back(offset);
        }
        in.close();
        set.FillBitmapWithBlockOffsets(block_offsets);
        set.Refresh_Block_Offsets();
        return true;
    }

    template<class T1, class T_FIELD, class SET_T>
    bool Load_Data(T_FIELD T1::* field, SET_T &set, const std::string &filename) {
        typedef  typename EnableIfSame<typename Allocator_type::template Array<T_FIELD>::type, T1, T>::type  Array_type;
        ifstream in(filename.c_str(), ios::in | ios::binary);
        if (!in.is_open())
            return false;
        set.Refresh_Block_Offsets();
        const unsigned long *const blocks = set.Get_Blocks().first;
        const unsigned numBlocks = set.Get_Blocks().second;
        const unsigned elementsPerBlock = alloc.Elements_Per_Block();
        Array_type data = alloc.Get_Array(field);
        T_FIELD *rawData = reinterpret_cast<T_FIELD*>(data.Get_Data_Ptr());
        for (int c = 0; c < numBlocks; c++) {
            in.read(reinterpret_cast<char*>((unsigned long)rawData + blocks[c]), elementsPerBlock * sizeof(T_FIELD));
        }
        in.close();
        return true;
    }

};
}
#endif
