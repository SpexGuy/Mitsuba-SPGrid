//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGrid_Set_h__
#define __SPGrid_Set_h__

#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <vector>
#include <pthread.h>

#include "std_array.h"

namespace SPGrid{

//#####################################################################
// Class SPGrid_Set
//#####################################################################
template<class T_ARRAY>
class SPGrid_Set
{
    enum {dim=T_ARRAY::dim};
    typedef typename T_ARRAY::DATA T;

public:
    typedef typename T_ARRAY::MASK T_MASK;
    T_ARRAY array;
    unsigned long* page_mask_array;
    unsigned long array_length;
    unsigned long max_linear_offset; // TODO: Change semantics to make this the first offset that is *invalid*
    std::vector<unsigned long> block_offsets;
    bool dirty; // Indicates that block offsets are inconsistent with the bitmap (perhaps for a good reason, if only one of them is used)

    // Single size argument constructor    
    SPGrid_Set(T_ARRAY array_input)
        :array(array_input)
    {
        unsigned long number_of_elements = array.geometry.Padded_Volume();
        unsigned long number_of_pages = number_of_elements>>T_MASK::block_bits;
        array_length = (number_of_pages+0x3fUL)>>6;
        page_mask_array=new unsigned long[array_length]; // TODO: Why "new" and not Raw_Allocate() ?
        memset(reinterpret_cast<void*>(page_mask_array),0,array_length*sizeof(unsigned long)); // TODO: Is this really needed?
        max_linear_offset = array.geometry.Padded_Volume()*T_ARRAY::MASK::Bytes_Per_Element(); // TODO: Check this is correct
        dirty = false; // They start as consistent -- both empty
    }

    SPGrid_Set(SPGrid_Set &&other)
        : array(std::move(other.array)),
          page_mask_array(std::move(other.page_mask_array)),
          array_length(std::move(other.array_length)),
          max_linear_offset(std::move(other.max_linear_offset)),
          block_offsets(std::move(other.block_offsets)),
          dirty(std::move(other.dirty))
    {}
    
    inline bool CheckBounds(unsigned long linear_offset)
    {
        return (linear_offset < max_linear_offset);
    }

    void MarkPageActive(unsigned long linear_offset)
    {
        if( linear_offset < max_linear_offset )
        {
            unsigned long page_mask = 1UL << (linear_offset>>12 & 0x3f);
            page_mask_array[linear_offset>>18] |= page_mask;
            dirty = true;
        } 
        else 
            FATAL_ERROR("Linear offset "+Value_To_String(linear_offset)+" is out of range (upper limit = "+Value_To_String(max_linear_offset)+")");
    }

    void Mask(const std_array<unsigned int,dim> coord,const T mask)
    {
        T& data=array(coord);
        unsigned long linear_offset=reinterpret_cast<unsigned long>(&data)-reinterpret_cast<unsigned long>(array.Get_Data_Ptr());
        unsigned long page_mask = 1UL << (linear_offset>>12 & 0x3f);
        assert((linear_offset>>18) < array_length);
        page_mask_array[linear_offset>>18] |= page_mask; dirty = true;
        data |= mask;

        // TODO: Examine using the following syntax instead
        // unsigned long linear_offset=T_ARRAY::MASK::Linear_Offset(coord);
        // MarkPageActive(linear_offset);
        // array(linear_offset) |= mask;
    }

    void Mask(unsigned long linear_offset, const T mask)
    {
        if( linear_offset < max_linear_offset)
        {
            T* data=reinterpret_cast<T*>(reinterpret_cast<unsigned long>(array.Get_Data_Ptr()) + linear_offset);
            unsigned long page_mask = 1UL << (linear_offset>>12 & 0x3f);
            page_mask_array[linear_offset>>18] |= page_mask; dirty = true;
            *data |= mask;
        } 
    }

    bool MaskMatch(unsigned long linear_offset, const T mask) const
    {
        if( linear_offset < max_linear_offset)
        {
            T* data=reinterpret_cast<T*>(reinterpret_cast<unsigned long>(array.Get_Data_Ptr()) + linear_offset);
            return (*data & mask);
        } 
        return false; // TODO: Do we need a FATAL ERROR?
    }

    bool Is_Set(unsigned long linear_offset,const T mask) const
    {
        if(linear_offset < max_linear_offset)
        {
            if(page_mask_array[linear_offset>>18] & (1UL<<(linear_offset>>12 & 0x3f)))
            {
                T* data=reinterpret_cast<T*>(reinterpret_cast<unsigned long>(array.Get_Data_Ptr()) + linear_offset);
                return (*data & mask);
            }
        }
        return false;        
    }
    
    bool Is_Set(std_array<int,dim> coord,const T mask) const
    {
        unsigned long linear_offset = T_MASK::Linear_Offset(coord);
        if(linear_offset < max_linear_offset)
        {
            if(page_mask_array[linear_offset>>18] & (1UL<<(linear_offset>>12 & 0x3f)))
            {
                T* data=reinterpret_cast<T*>(reinterpret_cast<unsigned long>(array.Get_Data_Ptr()) + linear_offset);
                return (*data & mask);
            }
        }
        return false;        
    }

    std::pair<const unsigned long*,unsigned> Get_Blocks() const
    {
        if(block_offsets.size())
            return std::pair<const unsigned long*,unsigned>(&block_offsets[0],block_offsets.size());
        else
            return std::pair<const unsigned long*,unsigned>(reinterpret_cast<const unsigned long *>(NULL),0);
    }

    void Refresh_Block_Offsets()
    {
        if(dirty)
            block_offsets = GenerateBlockOffsets();
        dirty = false;
    }

    std::vector<unsigned long> GenerateBlockOffsets()
    {
        std::vector<unsigned long> block_offsets;
        for (unsigned long i = 0; i<array_length; i++)
        {
            if(page_mask_array[i])
            {
                for (unsigned long pos=0; pos<64; pos++)
                {
                    if(page_mask_array[i] & (1UL<<pos))
                        block_offsets.push_back((i<<18)|(pos<<12));
                }
            }
        }
        return block_offsets;
    }

    void Clear_Bitmap()
    {
        memset(reinterpret_cast<void*>(page_mask_array),0,array_length*sizeof(unsigned long));
        dirty=true;
    }

    void Clear_Blocks()
    {
        std::vector<unsigned long>().swap(block_offsets);
        dirty=true;
    }

    void Clear()
    {Clear_Bitmap();Clear_Blocks();dirty=false;}

    void FillBitmapWithBlockOffsets(std::vector<unsigned long> block_offsets)
    {
        dirty = true;
        for (std::vector<unsigned long>::iterator it = block_offsets.begin(); it != block_offsets.end(); ++it)
        {
            unsigned long cur_offset = *it;
            page_mask_array[cur_offset>>18] |= ( 1UL << (cur_offset>>12 & 0x3f) );
        }
    }
};
}
#endif
