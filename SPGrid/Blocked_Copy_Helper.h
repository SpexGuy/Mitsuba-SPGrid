//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __Blocked_Copy_Helper__
#define __Blocked_Copy_Helper__
namespace SPGrid{

template<class T,int elements_per_block>
class Blocked_Copy_Helper
{
    T* const x;         // output stream
    const T* const y1;  // first input stream
    const T c;
    const T* const y2;  // second input stream
    const unsigned* const flags;  // flag stream
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    static const int prefetch_degree = 0;

public:
    explicit Blocked_Copy_Helper(T* const x_input,const T* const y_input, const T c_input,const T* const y2_input,const unsigned* const flags_input,const unsigned long* const b_input, const int size_input)
        :x(x_input),y1(y_input),c(c_input),y2(y2_input),flags(flags_input),b(b_input),size(size_input)
    {
    }

    void Run()
    {
        Run_Index_Range(0,size-1);
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};

template<class T,int elements_per_block>
class Blocked_Copy2_Helper
{
    T* const x;         // output stream
    const T* const y1;  // first input stream
    const T c;
    const unsigned* const flags;  // flag stream
    const unsigned long* const b;   // block offset stream
    const int size;     // number of blocks to process
    static const int prefetch_degree = 0;

public:
    explicit Blocked_Copy2_Helper(T* const x_input,const T* const y_input, const T c_input,const unsigned* const flags_input,const unsigned long* const b_input, const int size_input)
        :x(x_input),y1(y_input),c(c_input),flags(flags_input),b(b_input),size(size_input)
    {
    }

    void Run()
    {
        Run_Index_Range(0,size-1);
    }

//#####################################################################
    void Run_Parallel(const int number_of_partitions);
    void Run_Index_Range(const int index_start,const int index_end);
//#####################################################################
};
}
#endif
