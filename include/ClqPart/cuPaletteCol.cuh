#pragma once

#include <stdint.h>
#include <ClqPart/Types.h>
#include <ClqPart/graph.h>

void buildConfGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        NODE_T *d_confOffsets,
        NODE_T *d_confAdjList,
        NODE_T *d_nConflicts);

void buildCompGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        NODE_T *d_confOffsets,
        NODE_T *d_confAdjList,
        NODE_T *d_nConflicts);

template <typename OffsetTy>
void buildCooConfGraphDevice(const uint32_t *, const int, const NODE_T *, const NODE_T, const NODE_T, OffsetTy *, NODE_T *, OffsetTy *);

template <typename OffsetTy>
void buildCooConfGraphDevice(const uint32_t *, const int, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, OffsetTy *, NODE_T *, OffsetTy *);

template <typename OffsetTy>
void buildCooCompGraphDevice(const uint32_t *, const int, const NODE_T *, const NODE_T, const NODE_T, OffsetTy *, OffsetTy *);

template <typename OffsetTy>
void buildCsrConfGraphDevice(const NODE_T, const OffsetTy *, OffsetTy *, const NODE_T *, NODE_T *, const OffsetTy);

template <typename OffsetTy>
OffsetTy *cubExclusiveSum(const NODE_T, OffsetTy *);

// EGR / general-graph variant. Input is the conflict graph itself in CSR
// (d_inIA / d_inJA), so we only need to walk its existing edges and gate
// by colour-list overlap. Kernel uses one thread per directed edge with
// binary search of d_inIA to recover the source vertex; output is COO
// (Edge[]) plus per-vertex degree counts in d_confOffsets, exactly the
// shape that buildCsrConfGraphDevice expects downstream.
template <typename OffsetTy>
void buildEgrCooConfGraphDevice(
        const NODE_T *d_inIA,
        const NODE_T *d_inJA,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        const NODE_T num_dir_edges,
        OffsetTy *d_confOffsets,
        NODE_T *d_confAdjList,
        OffsetTy *d_nConflicts);

// Subset overload for recursive refinement. d_nodeListPos[v] is v's
// index in the active nodeList, or -1 if inactive. At recursion levels
// >= 1 d_colList is packed |nodeList| * n_colors, so the kernel
// indexes via the position, not the vertex ID.
template <typename OffsetTy>
void buildEgrCooConfGraphDevice(
        const NODE_T *d_inIA,
        const NODE_T *d_inJA,
        const NODE_T *d_nodeListPos,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        const NODE_T num_dir_edges,
        OffsetTy *d_confOffsets,
        NODE_T *d_confAdjList,
        OffsetTy *d_nConflicts);