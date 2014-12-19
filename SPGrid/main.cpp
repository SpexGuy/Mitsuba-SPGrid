//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "SPGrid_Allocator.h"
#include "SPGrid_Set.h"
#include "Blocked_Copy_Helper.h"
#include "Laplace_Helper.h"
#include "Save_Helper.h"
#include "Load_Helper.h"
#include "std_array.h"

#include "PTHREAD_QUEUE.h"
#include "ADAPTIVE_SPHERE_RASTERIZER.h"
#include "DENSE_CUBE_RASTERIZER.h"
#include "FACE_INITIALIZER.h"
#include "GEOMETRY_BLOCK.h"
#include "HIERARCHICAL_RASTERIZER.h"

//#define LOOP_AT_END
#define BLOCKED_COPY
#define DENSE_CUBE
#define LOAD_FLAGS
//#define LOAD_SMOKE
//#define SAVE_GRID

extern PTHREAD_QUEUE* pthread_queue;
using namespace SPGrid;

typedef float T;
typedef struct Foo_struct {
    T x,y,z;
    unsigned flags;
} Foo;
typedef SPGrid_Allocator<Foo,3> Foo_Allocator;
typedef SPGrid_Allocator<Foo,3>::Array<>::mask Foo_Mask;
typedef SPGrid_Allocator<Foo,3>::Array<T>::type Data_array_type;
typedef SPGrid_Allocator<Foo,3>::Array<const T>::type Const_data_array_type;
typedef SPGrid_Allocator<Foo,3>::Array<unsigned>::type Flags_array_type;
typedef SPGrid_Set<Flags_array_type> Flags_set_type;
typedef Save_Helper<Foo,3> Saver;
typedef Load_Helper<Foo,3> Loader;
typedef std_array<int,3> Vec3i;
typedef std_array<float,3> Vec3f;

int main(int argc,char* argv[]) {

    if (argc != 3) {
        printf("Please specify size (power of two), and number of threads\n");
        exit(1);
    }
    int size = atoi(argv[1]);
    if ((size & (size-1)) != 0) {
        printf("For this limited demo, size must be a power of two.\n");
        exit(1);
    }
    int n_threads = atoi(argv[2]);
    pthread_queue = new PTHREAD_QUEUE(n_threads);

#ifndef LOAD_SMOKE
#ifndef LOAD_FLAGS
    Foo_Allocator allocator(size,size,size);
    Data_array_type d1 = allocator.Get_Array(&Foo::x);
    Const_data_array_type d2 = allocator.Get_Const_Array(&Foo::y);
    Const_data_array_type d3 = allocator.Get_Const_Array(&Foo::z);
    Flags_array_type flags = allocator.Get_Array(&Foo::flags);
    Flags_set_type flag_set(flags);

    Vec3i imin(0);
    Vec3i imax(size);
    Vec3f Xmin(0.f);
    Vec3f Xmax(1.f);
    Vec3f center(.5f,.5f,.5f);
    float inner_radius=.3f;
    float outer_radius=.31f;

    int active_cells = 0;
#ifdef DENSE_CUBE
    if (size > 1024) {
        std::cout << "DENSE_CUBE mode has no memory savings, a size of " << size
                  << " will allocate roughly " 
                  << ((size>>10)*(size>>10)*(size>>10)*sizeof(Foo))
                  << "GB\nFeel free to remove this warning if needed.\n";
        exit(1);
    }
    std::cout << "Flagging active cells (in a dense cube)...";
    active_cells = DENSE_CUBE_RASTERIZER<Flags_set_type>::Rasterize(flag_set, size);
#else
    if (size < 256) {
        std::cout << "This is a sparse configuration, crank that size up!\n";
    }
    std::cout << "Flagging active cells (on narrow band)...";
    GEOMETRY_BLOCK block(imin,imax,Xmin,Xmax);
    ADAPTIVE_SPHERE_RASTERIZER<Flags_set_type> adaptive_sphere(flag_set,center,inner_radius,outer_radius);
    HIERARCHICAL_RASTERIZER<ADAPTIVE_SPHERE_RASTERIZER<Flags_set_type> > rasterizer(adaptive_sphere);
    rasterizer.Iterate(block);
    active_cells = adaptive_sphere.total_active;
#endif //DENSE_CUBE
    std::cout << "done.\n";
    uint64_t bigsize = size;
    std::cout << "Activated " << active_cells << " cells, out of a possible "
              << bigsize*bigsize*bigsize << "\n\n";

    flag_set.Refresh_Block_Offsets();
    // Face flag initialization
    FACE_INITIALIZER<Foo, 3>::Flag_Active_Faces(flag_set);
    printf("Finished flagging active cell faces.\n");

    T c = 0.5f;

#ifdef BLOCKED_COPY
    // Perform parallel x = y + (c * z)
    Blocked_Copy_Helper<T, Data_array_type::MASK::elements_per_block> helper(
        (T*)d1.Get_Data_Ptr(),
        (T*)d2.Get_Data_Ptr(),
        c,
        (T*)d3.Get_Data_Ptr(),
        (unsigned*)flags.Get_Data_Ptr(),
        flag_set.Get_Blocks().first,
        flag_set.Get_Blocks().second);
    helper.Run_Parallel(n_threads);
    printf("Finished running SAXPY kernel.\n");
#else // BLOCKED_COPY
    Laplace_Helper<T,NextLogTwo<sizeof(Foo)>::value,3> helper(
        (T*)d1.Get_Data_Ptr(),
        (T*)d3.Get_Data_Ptr(),
        (unsigned*)flags.Get_Data_Ptr(),
        flag_set.Get_Blocks().first,
        flag_set.Get_Blocks().second,
        1,
        2./3.);
    helper.Run_Parallel(n_threads);
    printf("Finished running Laplace kernel.\n");
#endif // BLOCKED_COPY

#else // LOAD_FLAGS
    Foo_Allocator temp = Loader::Load_Allocator("blocks.spmask");
    Foo_Allocator allocator(std::move(temp));
    Data_array_type d1 = allocator.Get_Array(&Foo::x);
    Const_data_array_type d2 = allocator.Get_Const_Array(&Foo::y);
    Const_data_array_type d3 = allocator.Get_Const_Array(&Foo::z);
    Flags_array_type flags = allocator.Get_Array(&Foo::flags);
    Flags_set_type flag_set(flags);

    Loader::Load_Mask(allocator, flag_set, "blocks.spmask");
    printf("Loaded mask from blocks.spmask\n");
    Loader::Load_Data(allocator, &Foo::flags, flag_set, "flags.spdata");
    printf("Loaded flags from flags.spdata\n");
    Loader::Load_Data(allocator, &Foo::x, flag_set, "density.spdata");
    printf("Loaded data from density.spdata\n");
#endif // LOAD_FLAGS

#else // LOAD_SMOKE
    ifstream input("smoke.vol", ios::binary | ios::in);
    input.seekg(0x18);
    const int xres = 128;
    const int yres = 128;
    const int zres = 50;
    
    float AABB[6];
    Static_Assert(sizeof(float) == 4);
    input.read(reinterpret_cast<char*>(AABB), 6*sizeof(float)); // read the AABB
    printf("AABB: (%f, %f, %f), (%f, %f, %f)\n", AABB[0], AABB[1], AABB[2], AABB[3], AABB[4], AABB[5]);
    
    Foo_Allocator allocator(128, 128, 64);
    Data_array_type d1 = allocator.Get_Array(&Foo::x);
    Const_data_array_type d2 = allocator.Get_Const_Array(&Foo::y);
    Const_data_array_type d3 = allocator.Get_Const_Array(&Foo::z);
    Flags_array_type flags = allocator.Get_Array(&Foo::flags);
    Flags_set_type flag_set(flags);
    int count = 0;
    int active = 0;
    float value;
    std_array<unsigned, 3> coord;
    for (coord(2) = 0; coord(2) < zres; coord(2)++) {
        for (coord(1) = 0; coord(1) < yres; coord(1)++) {
            for (coord(0) = 0; coord(0) < xres; coord(0)++) {
                if (!input.read(reinterpret_cast<char*>(&value), sizeof(float))) {
                    printf("FAIL AT (%d, %d, %d)!!\n", coord(0), coord(1), coord(2));
                    exit(-1);
                }
                if (value) {
                    flag_set.Mask(coord, 1U);
                    d1(coord) = value;
                    ++active;
                }
                ++count;
            }
        }
    }
    printf("Activated %d of %d cells\n", active, count);
#endif

#ifdef SAVE_GRID
    Saver::Save_Mask(allocator, flag_set, "blocks.spmask");
    Saver::Save_Data(allocator, &Foo::x, flag_set, "density.spdata");
    Saver::Save_Data(allocator, &Foo::flags, flag_set, "flags.spdata");
    printf("Finished saving.\n");
#endif
       
#ifdef LOOP_AT_END
    std::cout << "Looping forever, check out my memory usage!\n";
    while(1) {
        usleep(1000);
    }
#endif
} 

