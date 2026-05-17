#include <cuda.h>
#include <ClqPart/cuPaletteCol.cuh>
#include <cub/device/device_scan.cuh>


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
__global__ void build_conflict_graph_kernel(
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

__global__ void build_complement_graph_kernel(
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
                atomicAdd(d_nConflicts, 1);
                int index_offset = atomicAdd(&d_confOffsets[row], 1);
                d_confAdjList[row * n_vertices + index_offset] = col;
            }
        }
    }
}

template <typename OffsetTy>
__global__ void build_coo_complement_graph_kernel(
        const uint32_t *__restrict__ d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *__restrict__ d_colList, 
        const NODE_T n_vertices, 
        const NODE_T n_colors,
        OffsetTy *__restrict__ d_confOffsets,  
        OffsetTy *__restrict__ d_nConflicts){
    // int num_edges = n_vertices*(n_vertices - 1)/2;
    OffsetTy num_edges = (OffsetTy)n_vertices*(OffsetTy)n_vertices;
    // NODE_T halfway_point = n_vertices / 2;
    // Grid-Stride Loop only for lower triangle of matrix
    OffsetTy conflict_count = 0;
    for(OffsetTy edge_id = blockIdx.x * blockDim.x + threadIdx.x; edge_id < num_edges; edge_id += blockDim.x * gridDim.x){
        NODE_T row = edge_id / (OffsetTy)n_vertices;
        NODE_T col = edge_id - ((OffsetTy)row * (OffsetTy)n_vertices);
        if(row > col){
            const uint32_t *pauli1 = &d_pauliEnc[row * pauliEncSize];
            const uint32_t *pauli2 = &d_pauliEnc[col * pauliEncSize];
            bool isedge = compare_pauli_matrices(pauli1, pauli2, pauliEncSize);
            // If conflicting complement edge
            if(!isedge){
                conflict_count++;
                atomicAdd(&d_confOffsets[row], 1);
                atomicAdd(&d_confOffsets[col], 1);
            }
        }
    }
    atomicAdd(d_nConflicts, conflict_count);
}

template <typename OffsetTy>
__global__ void build_coo_conflict_graph_kernel(
        const uint32_t *__restrict__ d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *__restrict__ d_colList, 
        const NODE_T n_vertices, 
        const NODE_T n_colors,
        OffsetTy *__restrict__ d_confOffsets, 
        NODE_T *__restrict__ d_confAdjList, 
        OffsetTy *__restrict__ d_nConflicts){
    Edge *d_cooEdgeList = (Edge *)d_confAdjList;
    // int num_edges = n_vertices*(n_vertices - 1)/2;
    OffsetTy num_edges = (OffsetTy)n_vertices*(OffsetTy)n_vertices;
    // NODE_T halfway_point = n_vertices / 2;
    // Grid-Stride Loop only for lower triangle of matrix
    for(OffsetTy edge_id = blockIdx.x * blockDim.x + threadIdx.x; edge_id < num_edges; edge_id += blockDim.x * gridDim.x){
        NODE_T row = edge_id / (OffsetTy)n_vertices;
        NODE_T col = edge_id - ((OffsetTy)row * (OffsetTy)n_vertices);
        // Equivalent to edge_id % n_vertices
        if(row > col){
            const uint32_t *pauli1 = &d_pauliEnc[row * pauliEncSize];
            const uint32_t *pauli2 = &d_pauliEnc[col * pauliEncSize];
            bool isedge = compare_pauli_matrices(pauli1, pauli2, pauliEncSize);
            // If conflicting complement edge
            if(!isedge){
                const NODE_T *colList1 = &d_colList[row * n_colors];
                const NODE_T *colList2 = &d_colList[col * n_colors];
                bool common_color = findFirstCommonElement(colList1, colList2, n_colors);
                if(common_color){
                    OffsetTy index_offset = atomicAdd(d_nConflicts, 1);
                    atomicAdd(&d_confOffsets[row], 1);
                    atomicAdd(&d_confOffsets[col], 1);
                    d_cooEdgeList[index_offset] = Edge{row, col};
                }
            }
        }
    }
}

template <typename OffsetTy>
__global__ void build_coo_conflict_graph_kernel(
        const uint32_t *__restrict__ d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *__restrict__ d_colList, 
        const NODE_T *__restrict__ d_nodeList,
        const NODE_T n_vertices, 
        const NODE_T n_colors,
        OffsetTy *__restrict__ d_confOffsets, 
        NODE_T *__restrict__ d_confAdjList, 
        OffsetTy *__restrict__ d_nConflicts){
    Edge *d_cooEdgeList = (Edge *)d_confAdjList;
    // int num_edges = n_vertices*(n_vertices - 1)/2;
    OffsetTy num_edges = (OffsetTy)n_vertices*(OffsetTy)n_vertices;
    // NODE_T halfway_point = n_vertices / 2;
    // Grid-Stride Loop only for lower triangle of matrix
    for(OffsetTy edge_id = blockIdx.x * blockDim.x + threadIdx.x; edge_id < num_edges; edge_id += blockDim.x * gridDim.x){
        NODE_T row = edge_id / (OffsetTy)n_vertices;
        NODE_T col = edge_id - ((OffsetTy)row * (OffsetTy)n_vertices);
        // Equivalent to edge_id % n_vertices
        if(row > col){
            const NODE_T row_mapped = d_nodeList[row];
            const NODE_T col_mapped = d_nodeList[col];
            const uint32_t *pauli1 = &d_pauliEnc[row_mapped * pauliEncSize];
            const uint32_t *pauli2 = &d_pauliEnc[col_mapped * pauliEncSize];
            bool isedge = compare_pauli_matrices(pauli1, pauli2, pauliEncSize);
            // If conflicting complement edge
            if(!isedge){
                const NODE_T *colList1 = &d_colList[row * n_colors];
                const NODE_T *colList2 = &d_colList[col * n_colors];
                bool common_color = findFirstCommonElement(colList1, colList2, n_colors);
                if(common_color){
                    OffsetTy index_offset = atomicAdd(d_nConflicts, 1);
                    atomicAdd(&d_confOffsets[row_mapped], 1);
                    atomicAdd(&d_confOffsets[col_mapped], 1);
                    d_cooEdgeList[index_offset] = Edge{row_mapped, col_mapped};
                }
            }
        }
    }
}

template <typename OffsetTy>
__global__ void build_csr_conflict_graph_kernel(
        const NODE_T n_vertices, 
        const OffsetTy num_conf_edges,
        const OffsetTy *__restrict__ d_confOffsets, 
        OffsetTy *__restrict__ d_confOffsetsCnt, 
        const NODE_T *__restrict__ d_confAdjList,
        NODE_T *__restrict__ d_confCsr){
    Edge *d_cooEdgeList = (Edge *)d_confAdjList;
    for(OffsetTy edge_id = blockIdx.x * blockDim.x + threadIdx.x; edge_id < num_conf_edges; edge_id += blockDim.x * gridDim.x){
        Edge edge = d_cooEdgeList[edge_id];
        NODE_T u = edge.u;
        NODE_T v = edge.v;
        OffsetTy u_offset = atomicAdd(&d_confOffsetsCnt[u], 1);
        d_confCsr[d_confOffsets[u] + u_offset] = v;
        OffsetTy v_offset = atomicAdd(&d_confOffsetsCnt[v], 1);
        d_confCsr[d_confOffsets[v] + v_offset] = u;
    }
}

void buildCompGraphDevice(
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
    build_complement_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, n_vertices, n_colors, d_confOffsets, d_confAdjList, d_nConflicts);
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
    build_conflict_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, n_vertices, n_colors, d_confOffsets, d_confAdjList, d_nConflicts);
}

template <typename OffsetTy>
void buildCooConfGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        OffsetTy *d_confOffsets,
        NODE_T *d_confAdjList,
        OffsetTy *d_nConflicts){
    // Find cuda properties
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_coo_conflict_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, n_vertices, n_colors, d_confOffsets, d_confAdjList, d_nConflicts);
}

template <typename OffsetTy>
void buildCooConfGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T *d_nodeList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        OffsetTy *d_confOffsets,
        NODE_T *d_confAdjList,
        OffsetTy *d_nConflicts){
    // Find cuda properties
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_coo_conflict_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, d_nodeList, n_vertices, n_colors, d_confOffsets, d_confAdjList, d_nConflicts);
}

template <typename OffsetTy>
void buildCooCompGraphDevice(
        const uint32_t *d_pauliEnc,
        const int pauliEncSize,
        const NODE_T *d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        OffsetTy *d_confOffsets,
        OffsetTy *d_nConflicts){
    // Find cuda properties
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_coo_complement_graph_kernel<<<num_blocks, block_size>>>(d_pauliEnc, pauliEncSize, d_colList, n_vertices, n_colors, d_confOffsets, d_nConflicts);
}

template <typename OffsetTy>
void buildCsrConfGraphDevice(
        const NODE_T n_vertices,
        const OffsetTy *d_confOffsets,
        OffsetTy *d_confOffsetsCnt,
        const NODE_T *d_confAdjList,
        NODE_T *d_confCsr,
        const OffsetTy nConflicts){
    cudaMemset(d_confOffsetsCnt, 0, n_vertices * sizeof(OffsetTy));
    cudaDeviceSynchronize();
    // Call kernel to generate the CSR
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_csr_conflict_graph_kernel<<<num_blocks, block_size>>>(n_vertices, nConflicts, d_confOffsets, d_confOffsetsCnt, d_confAdjList, d_confCsr);
}

template <typename OffsetTy>
__host__ OffsetTy *cubExclusiveSum(const NODE_T n, OffsetTy *d_confOffsets){
    int num_elements = n + 1;
    size_t num_bytes = num_elements * sizeof(OffsetTy);
    // Find out how much temporary storage is needed for cub
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;
    cub::DeviceScan::ExclusiveSum(d_temp_storage, temp_storage_bytes, d_confOffsets, d_confOffsets, num_elements);
    std::cout << "Temp Storage Bytes: " << temp_storage_bytes << std::endl;
    std::cout << "Current Bytes: " << num_bytes << std::endl;
    size_t largest = std::max(temp_storage_bytes, num_bytes);
    cudaMalloc(&d_temp_storage, largest);
    // Call cub to perform exclusive sum
    cub::DeviceScan::ExclusiveSum<OffsetTy *, OffsetTy *>(d_temp_storage, largest, d_confOffsets, d_confOffsets, num_elements);
    return (OffsetTy *) d_temp_storage;
}

// EGR / LightGraph kernel. The input is already a conflict graph in CSR,
// so we just walk its directed-edge list and gate by colour-list overlap.
// One thread per directed edge; binary-search d_inIA to recover the
// source vertex.
template <typename OffsetTy>
__global__ void build_egr_coo_conflict_graph_kernel(
        const NODE_T *__restrict__ d_inIA,
        const NODE_T *__restrict__ d_inJA,
        const NODE_T *__restrict__ d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        const NODE_T num_dir_edges,
        OffsetTy *__restrict__ d_confOffsets,
        NODE_T *__restrict__ d_confAdjList,
        OffsetTy *__restrict__ d_nConflicts){
    Edge *d_cooEdgeList = (Edge *)d_confAdjList;
    for(NODE_T edge_id = blockIdx.x * blockDim.x + threadIdx.x;
        edge_id < num_dir_edges;
        edge_id += blockDim.x * gridDim.x){
        // Binary search: smallest u such that d_inIA[u+1] > edge_id.
        NODE_T lo = 0;
        NODE_T hi = n_vertices;
        while(lo < hi){
            NODE_T mid = (lo + hi) >> 1;
            if(d_inIA[mid + 1] <= edge_id) lo = mid + 1;
            else                            hi = mid;
        }
        NODE_T u = lo;
        NODE_T v = d_inJA[edge_id];
        if(v <= u) continue; // each undirected edge once
        const NODE_T *colList1 = &d_colList[(size_t)u * n_colors];
        const NODE_T *colList2 = &d_colList[(size_t)v * n_colors];
        if(findFirstCommonElement(colList1, colList2, n_colors)){
            OffsetTy idx = atomicAdd(d_nConflicts, (OffsetTy)1);
            atomicAdd(&d_confOffsets[u], (OffsetTy)1);
            atomicAdd(&d_confOffsets[v], (OffsetTy)1);
            d_cooEdgeList[idx] = Edge{u, v};
        }
    }
}

// Subset variant. d_nodeListPos[v] gives v's position in the active nodeList
// (or -1 if v is not active). At recursion level >= 1, d_colList has the
// packed layout |nodeList| * n_colors (assignListColor(nodeList,offset)
// only writes the active vertices), so we have to look up colList via
// the position, not the vertex ID. d_nodeListPos doubles as the inSet
// filter: a -1 entry means "skip this edge".
template <typename OffsetTy>
__global__ void build_egr_coo_conflict_graph_kernel(
        const NODE_T *__restrict__ d_inIA,
        const NODE_T *__restrict__ d_inJA,
        const NODE_T *__restrict__ d_nodeListPos,
        const NODE_T *__restrict__ d_colList,
        const NODE_T n_vertices,
        const NODE_T n_colors,
        const NODE_T num_dir_edges,
        OffsetTy *__restrict__ d_confOffsets,
        NODE_T *__restrict__ d_confAdjList,
        OffsetTy *__restrict__ d_nConflicts){
    Edge *d_cooEdgeList = (Edge *)d_confAdjList;
    for(NODE_T edge_id = blockIdx.x * blockDim.x + threadIdx.x;
        edge_id < num_dir_edges;
        edge_id += blockDim.x * gridDim.x){
        NODE_T lo = 0;
        NODE_T hi = n_vertices;
        while(lo < hi){
            NODE_T mid = (lo + hi) >> 1;
            if(d_inIA[mid + 1] <= edge_id) lo = mid + 1;
            else                            hi = mid;
        }
        NODE_T u = lo;
        NODE_T v = d_inJA[edge_id];
        if(v <= u) continue;
        NODE_T uPos = d_nodeListPos[u];
        NODE_T vPos = d_nodeListPos[v];
        if(uPos < 0 || vPos < 0) continue;
        const NODE_T *colList1 = &d_colList[(size_t)uPos * n_colors];
        const NODE_T *colList2 = &d_colList[(size_t)vPos * n_colors];
        if(findFirstCommonElement(colList1, colList2, n_colors)){
            OffsetTy idx = atomicAdd(d_nConflicts, (OffsetTy)1);
            atomicAdd(&d_confOffsets[u], (OffsetTy)1);
            atomicAdd(&d_confOffsets[v], (OffsetTy)1);
            d_cooEdgeList[idx] = Edge{u, v};
        }
    }
}

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
        OffsetTy *d_nConflicts){
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_egr_coo_conflict_graph_kernel<<<num_blocks, block_size>>>(
        d_inIA, d_inJA, d_colList, n_vertices, n_colors, num_dir_edges,
        d_confOffsets, d_confAdjList, d_nConflicts);
}

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
        OffsetTy *d_nConflicts){
    int device;
    cudaDeviceProp prop;
    cudaGetDevice(&device);
    cudaGetDeviceProperties(&prop, device);
    int nSM = prop.multiProcessorCount;
    int maxThreadsPerSM = prop.maxThreadsPerMultiProcessor;
    int block_size = 256;
    int num_blocks = nSM * (maxThreadsPerSM / block_size);
    build_egr_coo_conflict_graph_kernel<<<num_blocks, block_size>>>(
        d_inIA, d_inJA, d_nodeListPos, d_colList, n_vertices, n_colors, num_dir_edges,
        d_confOffsets, d_confAdjList, d_nConflicts);
}

// Create forced instantiation of templates for NODE_T, unsigned int, and unsigned long long
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, unsigned int *, NODE_T *, unsigned int *);
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, unsigned long long *, NODE_T *, unsigned long long *);
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, NODE_T *, NODE_T *, NODE_T *);
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, unsigned int *, NODE_T *, unsigned int *);
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, unsigned long long *, NODE_T *, unsigned long long *);
template void buildCooConfGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, NODE_T *, NODE_T *, NODE_T *);
template void buildCooCompGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, unsigned int *, unsigned int *);
template void buildCooCompGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, unsigned long long *, unsigned long long *);
template void buildCooCompGraphDevice(const unsigned int *, const int, const NODE_T *, const NODE_T, const NODE_T, NODE_T *, NODE_T *);
template void buildCsrConfGraphDevice(const NODE_T, const unsigned int *, unsigned int *, const NODE_T *, NODE_T *, const unsigned int);
template void buildCsrConfGraphDevice(const NODE_T, const unsigned long long *, unsigned long long *, const NODE_T *, NODE_T *, const unsigned long long);
template void buildCsrConfGraphDevice(const NODE_T, const NODE_T *, NODE_T *, const NODE_T *, NODE_T *, const NODE_T);
template unsigned int * cubExclusiveSum(const NODE_T, unsigned int *);
template unsigned long long * cubExclusiveSum(const NODE_T, unsigned long long *);
template NODE_T * cubExclusiveSum(const NODE_T, NODE_T *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, unsigned int *, NODE_T *, unsigned int *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, unsigned long long *, NODE_T *, unsigned long long *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, NODE_T *, NODE_T *, NODE_T *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, unsigned int *, NODE_T *, unsigned int *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, unsigned long long *, NODE_T *, unsigned long long *);
template void buildEgrCooConfGraphDevice(const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T *, const NODE_T, const NODE_T, const NODE_T, NODE_T *, NODE_T *, NODE_T *);
