#ifndef __Save_Helper__
#define __Save_Helper__

#include <iostream>
#include <fstream>
#include "SPGrid_Utilities.h"

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

namespace SPGrid{

template<class ALLOC_T, unsigned int dim, class SET_T> // TODO: enforce SET_T comes from ALLOC_T, dim == alloc_t.dim
struct Save_Helper
{
public:
    const ALLOC_T &alloc;
    SET_T &set;

    Save_Helper(const ALLOC_T &alloc_in, SET_T &set_in)
        : alloc(alloc_in), set(set_in)
    {}

    bool Save_Mask(const std::string &filename) {
        ofstream out(filename.c_str(), ios::out | ios::trunc);
        if (!out.is_open())
            return false;
        // first line is the block size
        std_array<unsigned, dim> size = alloc.Block_Size();
        out << size(0);
        for (int c = 1; c < dim; c++) {
            out << ' ' << size(c);
        }
        out << endl;
        // then the offsets
        set.Refresh_Block_Offsets();
        const unsigned long *blocks = set.Get_Blocks().first;
        unsigned numBlocks = set.Get_Blocks().second;
        for (int c = 0; c < numBlocks; c++) {
            out << blocks[c] << endl;
        }
        out.close();
        return true;
    }

    template<class ARRAY_T> // TODO: somehow ensure that the set and array come from the same allocator type
    bool Save_Data(const ARRAY_T &data, const std::string &filename) {
        typedef typename ARRAY_T::DATA T;
        ofstream out(filename.c_str(), ios::out | ios::trunc | ios::binary);
        if (!out.is_open())
            return false;
        set.Refresh_Block_Offsets();
        const unsigned long *blocks = set.Get_Blocks().first;
        const unsigned numBlocks = set.Get_Blocks().second;
        const unsigned elementsPerBlock = alloc.Elements_Per_Block();
        const T *rawData = reinterpret_cast<const T*>(data.Get_Data_Ptr());
        for (int c = 0; c < numBlocks; c++) {
            out.write(reinterpret_cast<const char*>((unsigned long)rawData + blocks[c]), elementsPerBlock * sizeof(T));
        }
        out.close();
        return true;
    }

};
}
#endif
