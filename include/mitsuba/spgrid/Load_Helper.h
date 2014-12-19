#ifndef __Load_Helper__
#define __Load_Helper__

#include <exception>
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
    static Allocator_type Load_Allocator(const std::string &filename) {
        ifstream in(filename.c_str(), ios::in);
        if (!in.is_open()) {
            std::cout << "Couldn't open " << filename << std::endl;
            throw std::exception();
        }

        // first line is the block size - use it to make the vector
        std_array<unsigned, dim> size;
        LoadSize(in, size);
        Allocator_type alloc(size);

        return alloc;
    }

    template<class SET_T>
    static void Load_Mask(Allocator_type &alloc, SET_T &set, const std::string &filename) {
        ifstream in(filename.c_str(), ios::in);
        if (!in.is_open()) {
            std::cout << "Couldn't open " << filename << std::endl;
            throw std::exception();
        }
        // first line is the block size - verify that it is the same
        std_array<unsigned, dim> fileSize;
        std_array<unsigned, dim> size = alloc.Size();
        LoadSize(in, fileSize);
        if (size != fileSize) {
            std::cout << "Bad grid size - expected " << size << " but was " << fileSize << std::endl;
            throw std::exception();
        }

        // then the offsets
        LoadSet(in, set);
    }

    template<class T1, class T_FIELD, class SET_T>
    static void Load_Data(Allocator_type &alloc, T_FIELD T1::* field, SET_T &set, const std::string &filename) {
        typedef  typename EnableIfSame<typename Allocator_type::template Array<T_FIELD>::type, T1, T>::type  Array_type;
        ifstream in(filename.c_str(), ios::in | ios::binary);
        if (!in.is_open()) {
            std::cout << "Couldn't open " << filename << std::endl;
            throw std::exception();
        }
        set.Refresh_Block_Offsets();
        const unsigned long *const blocks = set.Get_Blocks().first;
        const unsigned numBlocks = set.Get_Blocks().second;
        const unsigned elementsPerBlock = alloc.Elements_Per_Block();
        Array_type data = alloc.Get_Array(field);
        T_FIELD *rawData = reinterpret_cast<T_FIELD*>(data.Get_Data_Ptr());
        for (unsigned c = 0; c < numBlocks; c++) {
            in.read(reinterpret_cast<char*>((unsigned long)rawData + blocks[c]), elementsPerBlock * sizeof(T_FIELD));
        }
    }

private:
    template<class SET_T>
    static void LoadSet(std::istream &in, SET_T &set) {
        unsigned long offset;
        std::vector<unsigned long> block_offsets;
        while (in >> offset) {
            block_offsets.push_back(offset);
        }
        set.FillBitmapWithBlockOffsets(block_offsets);
        set.Refresh_Block_Offsets();
    }

    static void LoadSize(std::istream &in, std_array<unsigned, dim> &size) {
        std::string dims;
        if (!std::getline(in, dims)) {
            std::cout << "Empty file!" << std::endl;
            throw std::exception();
        }
        std::istringstream firstLine(dims);
        unsigned length;
        unsigned dimCount = 0;
        while(firstLine >> length) {
            if (dimCount >= dim) {
                std::cout << "Too many dimensions!" << std::endl;
                throw std::exception();
            }
            size(dimCount) = length;
            ++dimCount;
        }
        if (dimCount != dim) {
            std::cout << "Wrong number of dimensions! Expected " << dim << " but was " << dimCount << std::endl;
            throw std::exception();
        }
    }

};
}
#endif
