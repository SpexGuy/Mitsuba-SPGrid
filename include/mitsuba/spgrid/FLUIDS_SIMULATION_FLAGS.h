//#####################################################################
// Copyright (c) 2014, the authors of submission papers_0203
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#ifndef __FLUIDS_SIMULATION_FLAGS_h__
#define __FLUIDS_SIMULATION_FLAGS_h__

enum {

    // Cell properties

    SPGrid_Cell_Type_Interior  = 0x00000001u, // Set by user (also means this carries density)
    SPGrid_Cell_Type_Dirichlet = 0x00000002u, // Set by user
    SPGrid_Cell_Type_Ghost     = 0x00000004u, // Generated automatically
    SPGrid_Cell_Type_Active    = 0x00000008u, // Active = Interior but not Dirichlet, also not fully surrounded by Neumann faces

    // Face properties

    SPGrid_Face_Type_X_Valid   = 0x00000010u, // Generated automatically - A face is valid if :
    SPGrid_Face_Type_Y_Valid   = 0x00000020u, //   (a) it is the face of an interior cell, AND
    SPGrid_Face_Type_Z_Valid   = 0x00000040u, //   (b) is undivided (not a part of a larger face)

    SPGrid_Face_Type_X_Active  = 0x00000080u, // Active faces are initialized to be all the valid faces that don't touch exterior cells
    SPGrid_Face_Type_Y_Active  = 0x00000100u, // These flags can be overridden by user, to specify additional Neumann regions
    SPGrid_Face_Type_Z_Active  = 0x00000200u, // NOTE: The user can only *deactivate* active faces, not vice versa

    // Node properties

    SPGrid_Node_Active         = 0x00000400u,
    SPGrid_Node_Coarse_Shared  = 0x00000800u,
    
    // Determined by cell types
    SPGrid_Face_Minus_X_Active = 0x00001000u,
    SPGrid_Face_Plus_X_Active  = 0x00002000u,
    SPGrid_Face_Minus_Y_Active = 0x00004000u,
    SPGrid_Face_Plus_Y_Active  = 0x00008000u,
    SPGrid_Face_Minus_Z_Active = 0x00010000u,
    SPGrid_Face_Plus_Z_Active  = 0x00020000u,

    SPGrid_Face_Minus_X_Scaled = 0x00040000u,
    SPGrid_Face_Plus_X_Scaled  = 0x00080000u,
    SPGrid_Face_Minus_Y_Scaled = 0x00100000u,
    SPGrid_Face_Plus_Y_Scaled  = 0x00200000u,
    SPGrid_Face_Minus_Z_Scaled = 0x00400000u,
    SPGrid_Face_Plus_Z_Scaled  = 0x00800000u,

    SPGrid_Ghost_Child_000     = 0x01000000u,
    SPGrid_Ghost_Child_001     = 0x02000000u,
    SPGrid_Ghost_Child_010     = 0x04000000u,
    SPGrid_Ghost_Child_011     = 0x08000000u,
    SPGrid_Ghost_Child_100     = 0x10000000u,
    SPGrid_Ghost_Child_101     = 0x20000000u,
    SPGrid_Ghost_Child_110     = 0x40000000u,
    SPGrid_Ghost_Child_111     = 0x80000000u
};
#endif
