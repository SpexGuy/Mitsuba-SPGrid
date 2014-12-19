//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __SPGrid_Allocator_h__
#define __SPGrid_Allocator_h__

#include "SPGrid_Array.h"
#include "SPGrid_Mask.h"
#include "SPGrid_Utilities.h"

namespace SPGrid{

//#####################################################################
// Class SPGrid_Allocator
//#####################################################################
template<class T,int dim>
class SPGrid_Allocator: public SPGrid_Geometry<dim>, public SPGrid_Mask_base<NextLogTwo<sizeof(T)>::value,dim>
{
    typedef SPGrid_Geometry<dim> T_Geometry_base;
    typedef SPGrid_Mask_base<NextLogTwo<sizeof(T)>::value,dim> T_Mask_base;
    using T_Mask_base::block_bits;using T_Mask_base::block_xbits;using T_Mask_base::block_ybits;

private:
    void* data_ptr;

public:
    using T_Geometry_base::Padded_Volume;
    
    template<class T_FIELD=T> struct Array
    {
      typedef SPGrid_Mask<NextLogTwo<sizeof(T)>::value,NextLogTwo<sizeof(T_FIELD)>::value,dim> mask;
      typedef SPGrid_Array<T_FIELD,mask> type;
    };

    SPGrid_Allocator(const unsigned xsize_input, const unsigned ysize_input, const unsigned zsize_input)
        :SPGrid_Geometry<dim>(xsize_input,ysize_input,zsize_input,block_xbits,block_ybits,T_Mask_base::block_zbits)
    {
        Static_Assert(dim==3);
        Check_Compliance();
        data_ptr=Raw_Allocate(Padded_Volume()*Next_Power_Of_Two(sizeof(T)));
    }

    SPGrid_Allocator(const unsigned xsize_input, const unsigned ysize_input)
        :SPGrid_Geometry<dim>(xsize_input,ysize_input,block_xbits,block_ybits)
    {
        Static_Assert(dim==2);
        Check_Compliance();
        data_ptr=Raw_Allocate(Padded_Volume()*Next_Power_Of_Two(sizeof(T)));
    }
    
    SPGrid_Allocator(const std_array<unsigned int,3> size_in)
        :SPGrid_Geometry<dim>(size_in.data[0],size_in.data[1],size_in.data[2],block_xbits,block_ybits,T_Mask_base::block_zbits)
    {
        Static_Assert(dim==3);
        Check_Compliance();
        data_ptr=Raw_Allocate(Padded_Volume()*Next_Power_Of_Two(sizeof(T)));
    }

    SPGrid_Allocator(const std_array<unsigned int,2> size_in)
        :SPGrid_Geometry<dim>(size_in.data[0],size_in.data[1],block_xbits,block_ybits)
    {
        Static_Assert(dim==2);
        Check_Compliance();
        data_ptr=Raw_Allocate(Padded_Volume()*Next_Power_Of_Two(sizeof(T)));
    }

    SPGrid_Allocator(SPGrid_Allocator &&that)
        :SPGrid_Geometry<dim>(std::move(that)),
         data_ptr(std::move(that.data_ptr))
    {
        that.data_ptr = NULL;
    }

    ~SPGrid_Allocator()
    {
        if (data_ptr) {
            printf("Deallocating Live SPGrid_Allocator\n");
            Raw_Deallocate(data_ptr,Padded_Volume()*Next_Power_Of_Two(sizeof(T)));
        } else {
            printf("Deallocating null SPGrid_Allocator\n");
        }
    }

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<T_FIELD>::type,T1,T>::type
    Get_Array(T_FIELD T1::* field)
    {
        typedef typename Array<T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<unsigned long>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<const T_FIELD>::type,T1,T>::type
    Get_Array(T_FIELD T1::* field) const
    {
        typedef typename Array<const T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<unsigned long>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<const T_FIELD>::type,T1,T>::type
    Get_Const_Array(T_FIELD T1::* field) const
    {
        typedef typename Array<const T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<unsigned long>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    typename Array<>::type Get_Array()
    {
        typedef typename Array<>::type Array_type;
        return Array_type(data_ptr,*this);
    }

    typename Array<const T>::type Get_Array() const
    {
        typedef typename Array<const T>::type Array_type;
        return Array_type(data_ptr,*this);
    }

    typename Array<const T>::type Get_Const_Array() const
    {
        typedef typename Array<const T>::type Array_type;
        return Array_type(data_ptr,*this);
    }

protected:    
    inline void* Get_Data_Ptr() {return data_ptr;}

//#####################################################################
};
}
#endif
