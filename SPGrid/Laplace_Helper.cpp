//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include "Laplace_Helper.h"
#include "FLUIDS_SIMULATION_FLAGS.h"

using namespace SPGrid;

template<class T,int log2_struct,int d>
struct Laplace_Thread_Helper:public PTHREAD_QUEUE::TASK
{
    Laplace_Helper<T,log2_struct,d>* const obj;
    const int index_start,index_end;
    Laplace_Thread_Helper(Laplace_Helper<T,log2_struct,d>* const obj_input,const int index_start_input,const int index_end_input)
        :obj(obj_input),index_start(index_start_input),index_end(index_end_input) {}
    void Run(){obj->Run_Index_Range(index_start,index_end);}
};

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Laplace_Helper<T,log2_struct,2>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Laplace_Thread_Helper<T,log2_struct,d>* task=new Laplace_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}

template <class T, int log2_struct>
void Laplace_Helper<T,log2_struct,2>::ComputeShadowGrid( unsigned long* offset_grid_ptr, 
                                                         unsigned* mask_in, 
                                                         unsigned long packed_offset)
{
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    unsigned long simple_offset = 0;
    // Fill in simple offsets
    for (int i = xmin; i<=xmax; i++)
    for (int j = ymin; j<=ymax; j++)
    {
      o_grid[i][j] = packed_offset + simple_offset;  // Can do simple addition here since addresses are within block
      simple_offset += sizeof(T);
    }

    // First let's do the starting points
    o_grid[xmin-1][ymin] =  T_MASK::template Packed_OffsetXdim<-1>(o_grid[xmin][ymin]);
    o_grid[xmax+1][ymin] =  T_MASK::template Packed_OffsetXdim< 1>(o_grid[xmax][ymin]);

    o_grid[xmin][ymin-1] =  T_MASK::template Packed_OffsetYdim<-1>(o_grid[xmin][ymin]);
    o_grid[xmin][ymax+1] =  T_MASK::template Packed_OffsetYdim< 1>(o_grid[xmin][ymax]);

    // Fill in edge offsets (cube faces, but not edges will be correct after this)
    // This is ok for 6 neighbors, but one more pass will be needed for kernels that use edges
    {
      // Left and Right face
      for (int i=xmin-1; i<=xmax+1; i+= (xmax-xmin)+2)
      {
        simple_offset = o_grid[i][ymin];
        for (int j=ymin; j<=ymax; j++)
        {
          o_grid[i][j] = simple_offset;
          simple_offset += sizeof(T);  // Simple addition (going through neighboring block in same manner)
        }
      }
    }
    
    {
      // Top and bottom face
      for (int j = ymin-1; j<=ymax+1; j+= (ymax-ymin)+2)
      {
        simple_offset = o_grid[xmin][j];
        for (int i=xmin; i<=xmax; i++)
        {
          o_grid[i][j] = simple_offset;
          simple_offset += sizeof(T) * (block_ysize);
        }
      }
    }

}
//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Laplace_Helper<T,log2_struct,2>::Run_Index_Range(const int index_start,const int index_end)
{  
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* y1_in  = reinterpret_cast<T*>((unsigned long)y + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)y1_in - (unsigned long)y;
      unsigned long y_base_addr = reinterpret_cast<unsigned long>(y);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)) {
            double result=(double)0.;

            if (mask_value & SPGrid_Face_Minus_X_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Minus_X_Active)
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j])) - y1_in[cur_index]);

            if (mask_value & SPGrid_Face_Plus_X_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Plus_X_Active)
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j])) - y1_in[cur_index]);

            
            if (mask_value & SPGrid_Face_Minus_Y_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Minus_Y_Active)                                
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1])) - y1_in[cur_index]);

            if (mask_value & SPGrid_Face_Plus_Y_Scaled)                                      
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Plus_Y_Active)                                 
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1])) - y1_in[cur_index]);
            
            output[cur_index] = (T)result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
}

//#####################################################################
// Function Run_Parallel
//#####################################################################
template<class T, int log2_struct> void Laplace_Helper<T,log2_struct,3>::Run_Parallel(const int number_of_partitions)
{
    for(int partition=0;partition<number_of_partitions;partition++){
        // Calculate indicies of current partition
        int first_index_of_partition=(size/number_of_partitions)*(partition)+std::min(size%number_of_partitions,partition);
        int last_index_of_partition=(size/number_of_partitions)*(partition+1)+std::min(size%number_of_partitions,partition+1)-1;
        // Create helper object
        Laplace_Thread_Helper<T,log2_struct,d>* task=new Laplace_Thread_Helper<T,log2_struct,d>(this,first_index_of_partition,last_index_of_partition);
        // Enqueue
        pthread_queue->Queue(task);
    }
    // Wait for all tasks to complete
    pthread_queue->Wait();
}


template <class T, int log2_struct>
void Laplace_Helper<T,log2_struct,3>::ComputeShadowGrid( unsigned long* offset_grid_ptr, 
                                                                            unsigned* mask_in, 
                                                                            unsigned long packed_offset)
{
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    unsigned long simple_offset = 0;
    // Fill in simple offsets
    for (int i = xmin; i<=xmax; i++)
    for (int j = ymin; j<=ymax; j++)
    for (int k = zmin; k<=zmax; k++)
    {
      o_grid[i][j][k] = packed_offset + simple_offset;  // Can do simple addition here since addresses are within block
      simple_offset += sizeof(T);
    }

    // First let's do the starting points
    o_grid[xmin-1][ymin][zmin] =  T_MASK::template Packed_OffsetXdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmax+1][ymin][zmin] =  T_MASK::template Packed_OffsetXdim< 1>(o_grid[xmax][ymin][zmin]);

    o_grid[xmin][ymin][zmin-1] =  T_MASK::template Packed_OffsetZdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmin][ymin][zmax+1] =  T_MASK::template Packed_OffsetZdim< 1>(o_grid[xmin][ymin][zmax]);

    o_grid[xmin][ymin-1][zmin] =  T_MASK::template Packed_OffsetYdim<-1>(o_grid[xmin][ymin][zmin]);
    o_grid[xmin][ymax+1][zmin] =  T_MASK::template Packed_OffsetYdim< 1>(o_grid[xmin][ymax][zmin]);

    // Fill in edge offsets (cube faces, but not edges will be correct after this)
    // This is ok for 6 neighbors, but one more pass will be needed for kernels that use edges
    {
      // Left and Right face
      for (int i=xmin-1; i<=xmax+1; i+= (xmax-xmin)+2)
      {
        simple_offset = o_grid[i][ymin][zmin];
        for (int j=ymin; j<=ymax; j++)
        for (int k=zmin; k<=zmax; k++)
        {
          o_grid[i][j][k] = simple_offset;
          simple_offset += sizeof(T);  // Simple addition (going through neighboring block in same manner)
        }
      }
    }

    {
      // Front and Back face
      for (int k=zmin-1; k<=zmax+1; k+= (zmax-zmin)+2)
      {
        simple_offset = o_grid[xmin][ymin][k];
        for (int i=xmin; i<=xmax; i++)
        for (int j=ymin; j<=ymax; j++)
        {
          o_grid[i][j][k] = simple_offset;
          simple_offset += block_zsize*sizeof(T);  
        }
      }
    }
    
    {
      // Top and bottom face
      for (int j=ymin-1; j<=ymax+1; j+= (ymax-ymin)+2)
      {
        simple_offset = o_grid[xmin][j][zmin];
        for (int i=xmin; i<=xmax; i++)
        {
          for (int k=zmin; k<=zmax; k++)
          {
            o_grid[i][j][k] = simple_offset;
            simple_offset += sizeof(T);  
          }
          simple_offset += sizeof(T) * (block_ysize-1) * (block_zsize);
        }
      }
    }

}

//#####################################################################
// Function Run_Index_Range
//#####################################################################
// T_MASK corresponds to the mask for the data (not the mask channel)
template <class T, int log2_struct> void Laplace_Helper<T,log2_struct,3>::Run_Index_Range(const int index_start,const int index_end)
{  
    int count = 0;
    // Compute shadow grid of linear offsets
    unsigned long* offset_grid_ptr = (unsigned long*)malloc( (og_xsize) * (og_ysize) * (og_zsize) * sizeof(unsigned long));
    typedef unsigned long (&offset_grid_type)[og_xsize][og_ysize][og_zsize];
    offset_grid_type o_grid = reinterpret_cast<offset_grid_type>(*offset_grid_ptr);
    
    // Iterating through all block indices
    for(int index=index_start;index<=index_end;index++) {
      T* output = reinterpret_cast<T*>((unsigned long)x + b[index]);
      T* y1_in  = reinterpret_cast<T*>((unsigned long)y + b[index]);

      unsigned* mask_in  = reinterpret_cast<unsigned*>((unsigned long)mask + b[index]);
      unsigned long packed_offset = (unsigned long)y1_in - (unsigned long)y;
      unsigned long y_base_addr = reinterpret_cast<unsigned long>(y);

      ComputeShadowGrid(offset_grid_ptr, mask_in, packed_offset);

      int cur_index = 0;
      // Actually process elements
      for(int i=xmin;i<=xmax;i++)
      for(int j=ymin;j<=ymax;j++)
      for(int k=zmin;k<=zmax;k++)
      {
        unsigned mask_value = mask_in[cur_index];

        if ( mask_value & (SPGrid_Cell_Type_Interior|SPGrid_Cell_Type_Ghost)) {
            count++;
            double result=(double)0.;

            if (mask_value & SPGrid_Face_Minus_X_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j][k])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Minus_X_Active)
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i-1][j][k])) - y1_in[cur_index]);

            if (mask_value & SPGrid_Face_Plus_X_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j][k])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Plus_X_Active)
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i+1][j][k])) - y1_in[cur_index]);

            
            if (mask_value & SPGrid_Face_Minus_Y_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1][k])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Minus_Y_Active)                                
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j-1][k])) - y1_in[cur_index]);

            if (mask_value & SPGrid_Face_Plus_Y_Scaled)                                      
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1][k])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Plus_Y_Active)                                 
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j+1][k])) - y1_in[cur_index]);
            
            
            if (mask_value & SPGrid_Face_Minus_Z_Scaled)
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k-1])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Minus_Z_Active)                                   
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k-1])) - y1_in[cur_index]);
            
            if (mask_value & SPGrid_Face_Plus_Z_Scaled)                                         
              result += scale_nonuniform * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k+1])) - y1_in[cur_index]);
            else if (mask_value & SPGrid_Face_Plus_Z_Active)                                    
              result += scale_uniform    * ((double)(*reinterpret_cast<T*>(y_base_addr + o_grid[i][j][k+1])) - y1_in[cur_index]);
            
            output[cur_index] = (T)result;
        }
        cur_index++;
      }
    }
    free(offset_grid_ptr);
    std::cout << "Laplace helper saw " << count << " active cells.";
}

template class Laplace_Helper<float,5,2>;
template class Laplace_Helper<float,6,2>;

template class Laplace_Helper<float,4,3>;
template class Laplace_Helper<float,5,3>;
template class Laplace_Helper<float,6,3>;
//template class Laplace_Helper<float,5,2>;
