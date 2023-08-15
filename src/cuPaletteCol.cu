#include "ClqPart/cuPaletteCol.cuh"
#include <cuda.h>

__device__ __inline__ bool compare_pauli_matrices(
        const uint32_t * __restrict__ pauli1,
        const uint32_t * __restrict__ pauli2,
        const int pauliSize){
    uint32_t cnt = 0;
    for(int i = 0; i < pauliSize; i++){
        cnt += __popc(pauli1[i] & pauli2[i]);
    }
    if (cnt & 0x1) {
      return true;
    }
    else {
      return false;
    }
 }

__device__ __inline__ bool findFirstCommonElement(
        const NODE_T * __restrict__ colList1,
        const NODE_T * __restrict__ colList2,
        const NODE_T colSize) {
    int i = 0; // Index for colList1
    int j = 0; // Index for colList2

    while (i < colSize && j < colSize) {
        if (colList1[i] < colList2[j]) {
            i++; // Move to the next element in colList1
        } else if (colList1[i] > colList2[j]) {
            j++; // Move to the next element in colList2
        } else {
            return true; // Found a common element
        }
    }

    return false; // No common element found
}


// extern __shared__ uint32_t shared[];
__global__ void build_conf_graph_kernel(
        const uint32_t *__restrict__ d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *__restrict__ d_colList, 
        const NODE_T n_vertices, 
        const NODE_T n_colors,
        NODE_T *__restrict__ d_confOffsets, 
        NODE_T *__restrict__ d_confAdjList, 
        NODE_T *__restrict__ d_nConflicts){
    // NODE_T *s_pauliEnc = (uint32_t *)shared;
    // NODE_T *s_colList = (NODE_T *)&s_pauliEnc[pauliEncSize * shared_edges_size];
    int num_edges = n_vertices*n_vertices;
    // int block_edges = shared_edges_size * shared_edges_size;
    // Grid-Stride Loop
    for(int i = blockIdx.x * blockDim.x + threadIdx.x; i < num_edges; i += blockDim.x * gridDim.x){
        int row = i / n_vertices;
        int col = i % n_vertices;
        if(row != col){
            const uint32_t *pauli1 = &d_pauliEnc[row * pauliEncSize];
            const uint32_t *pauli2 = &d_pauliEnc[col * pauliEncSize];
            bool isedge = compare_pauli_matrices(pauli1, pauli2, pauliEncSize);
            // If conflicting complement edge
            if(!isedge){
                const NODE_T *colList1 = &d_colList[row * n_colors];
                const NODE_T *colList2 = &d_colList[col * n_colors];
                bool common_color = findFirstCommonElement(colList1, colList2, n_colors);
                if(common_color){
                    atomicAdd(d_nConflicts, 1);
                    int index_offset = atomicAdd(&d_confOffsets[row], 1);
                    d_confAdjList[row * n_vertices + index_offset] = col;
                }
            }
        }
    }
}

void buildConfGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        NODE_T *d_confOffsets,
        NODE_T *d_confAdjList,
        NODE_T *d_nConflicts){
    // Find cuda properties
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_conf_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, n_vertices, n_colors, d_confOffsets, d_confAdjList, d_nConflicts);
}