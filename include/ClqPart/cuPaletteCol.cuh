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
void buildCooCompGraphDevice(const uint32_t *, const int, const NODE_T *, const NODE_T, const NODE_T, OffsetTy *, OffsetTy *);

template <typename OffsetTy>
void buildCsrConfGraphDevice(const NODE_T, const OffsetTy *, OffsetTy *, const NODE_T *, NODE_T *, const OffsetTy);

template <typename OffsetTy>
void cubInclusiveSum(void *, const NODE_T, OffsetTy *);