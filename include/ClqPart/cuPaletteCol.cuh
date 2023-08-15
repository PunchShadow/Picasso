#pragma once

#include <stdint.h>
#include <ClqPart/Types.h>

// __inline__ bool compare_pauli_matrices(
//         const uint32_t * __restrict__ pauli1,
//         const uint32_t * __restrict__ pauli2,
//         const int pauliSize);


// __inline__ bool findFirstCommonElement(
//         const NODE_T * __restrict__ colList1,
//         const NODE_T * __restrict__ colList2,
//         const NODE_T colSize);


// __global__ build_conf_graph_kernel(
//         const uint32_t *__restrict__ d_pauliEnc,
//         const int pauliEncSize,
//         const NODE_T *__restrict__ d_colList, 
//         const NODE_T n_vertices, 
//         const NODE_T n_colors,
//         NODE_T *__restrict__ d_confOffsets, 
//         NODE_T *__restrict__ d_confAdjList, 
//         NODE_T *__restrict__ d_nConflicts);


void buildConfGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        NODE_T *d_confOffsets,
        NODE_T *d_confAdjList,
        NODE_T *d_nConflicts);