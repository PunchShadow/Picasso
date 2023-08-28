/*
 * Copyright (C) 2023  Ferdous,S M <ferdous.csebuet@egmail.com>
 * Author: Ferdous,S M <ferdous.csebuet@egmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// #define ENABLE_GPU
// #define COMPLEMENT_GRAPH

#include "ClqPart/graph.h" 
#include "ClqPart/JsonGraph.h"
#include <random>

#include <omp.h>

#ifdef ENABLE_GPU
// #include <cuda.h>
#include <cuda_runtime.h>
#include <cstdio>
#include <stdlib.h>
#include <atomic>
// CUDA Error Checking Code
#define ERR_CHK(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line,
    bool abort = true) {
  if (code != cudaSuccess) {
    fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code),
        file, line);
    if (abort)
      exit(code);
  }
}
#include "ClqPart/cuPaletteCol.cuh"
#endif

bool findFirstCommonElement(const std::vector<NODE_T>& vec1, const std::vector<NODE_T>& vec2);

template<typename OffsetTy = int>
class PaletteColor {
  
  NODE_T n;
  NODE_T colThreshold;
  std::vector<std::vector<NODE_T> > colList;
  std::vector<NODE_T> colors;
  std::vector<NODE_T> confColors;
  NODE_T nColors;
  OffsetTy nConflicts;

  #ifdef ENABLE_GPU
  std::vector<NODE_T> h_colList;
  NODE_T *d_colList;
  NODE_T *d_colors;
  NODE_T *d_confColors;
  std::vector<uint32_t> h_pauliEnc;
  uint32_t *d_pauliEnc;
  size_t pauliEncSize;
  std::vector<NODE_T> h_confVertices;
  NODE_T *d_confVertices;
  std::vector<OffsetTy> h_confOffsets;
  OffsetTy *d_confOffsets;
  #endif
  
  std::vector<NODE_T> confVertices;
  std::vector<NODE_T> invalidVertices;
  std::vector<NODE_T> vertexOrder;
  std::vector<std::vector< NODE_T> > confAdjList;

  NODE_T T;
  double assignTime;
  double confColorTime;
  double invalidColorTime;


public:
  PaletteColor( NODE_T n1, NODE_T colThresh, float alpha=1) {
    n = n1;
    colThreshold = colThresh; 
    colors.resize(n,-1);
    confColors.resize(n,-1);
    vertexOrder.resize(n);
    nColors = 0;
    nConflicts = 0;
    T =  static_cast<NODE_T> (alpha*log(n));
    
    confVertices.resize(n,0);
    confAdjList.resize(n);
    colList.resize(n);
    
    invalidVertices.clear();
    assignListColor();
  }

  template<typename PauliTy = std::string>
  void naiveGreedyColor(std::vector<NODE_T> vertList, ClqPart::JsonGraph &jsongraph, NODE_T offset) {

    if(vertList.empty() == false) {

      std::vector<NODE_T> forbiddenCol(n,-1);
      colors[vertList[0]] = offset; 

      for(auto i=1; i<vertList.size();i++) {
        NODE_T eu = vertList[i]; 
        for(auto j=0; j<i; j++) {
          NODE_T ev = vertList[j]; 

          if(jsongraph.is_an_edge<PauliTy>(eu,ev) == false) { 
            if (colors[ev] >= 0) {
              forbiddenCol[colors[ev]] = eu;
            }
          } 
        }
      
        //color eu with first available color
        for( auto ii=offset;ii<n;ii++) {
          if(forbiddenCol[ii] == -1) {
            colors[eu] = ii; 
            break;
          } 
        }
      }
      nColors = *std::max_element(colors.begin(),colors.end()) + 1;
    }
  }

  //This function computes the graph directly, rather in streaming way. It takes
//a JsonGraph object since it requires to determine whether the pair (eu,ev) is 
//an edge in the complement graph.
template<typename PauliTy = std::string>
void buildConfGraph ( ClqPart::JsonGraph &jsongraph) {

  
  for(NODE_T eu =0; eu < n-1; eu++) {
    for(NODE_T ev = eu+1; ev < n; ev++) {

      if(jsongraph.is_an_edge<PauliTy>(eu,ev) == false) {
        bool hasCommon = findFirstCommonElement(colList[eu],colList[ev]);
        if(hasCommon == true ) {

          confAdjList[eu].push_back(ev); 
          confAdjList[ev].push_back(eu); 
          nConflicts++;
        }
      }
    }
  }

}


//This function checks whether the edge (eu,ev) (in complement graph) is a conflict. NOte that 
//this is for constructing the graph in streaming way, so we are presented with
//one edge at a time. (eu,ev) is an edge in a complement graph.
void buildStreamConfGraph ( NODE_T eu, NODE_T ev) {

  bool hasCommon = findFirstCommonElement(colList[eu],colList[ev]);
  if(hasCommon == true ) {

    confAdjList[eu].push_back(ev); 
    confAdjList[ev].push_back(eu); 
    nConflicts++;
  }

}

#ifdef ENABLE_GPU
void buildConfGraphGpuMemConscious (ClqPart::JsonGraph &jsongraph) {
  pauliEncSize = jsongraph.getPauliEncSize();
  std::vector<uint32_t> &h_pauliEnc = jsongraph.getEncodedData();
  ERR_CHK(cudaMalloc(&d_pauliEnc, h_pauliEnc.size() * sizeof(NODE_T)));
  ERR_CHK(cudaMemcpy(d_pauliEnc, h_pauliEnc.data(), h_pauliEnc.size() * sizeof(NODE_T), cudaMemcpyHostToDevice));

  h_confOffsets.resize(n+1);
  ERR_CHK(cudaMalloc(&d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy)));
  ERR_CHK(cudaMemset(d_confOffsets, 0, h_confOffsets.size() * sizeof(OffsetTy)));
  OffsetTy *d_confOffsetsCnt;
  ERR_CHK(cudaMalloc(&d_confOffsetsCnt, n * sizeof(OffsetTy)));

  OffsetTy *d_nConflicts;
  ERR_CHK(cudaMalloc(&d_nConflicts, sizeof(OffsetTy)));
  ERR_CHK(cudaMemset(d_nConflicts, 0, sizeof(OffsetTy)));

  #ifndef COMPLEMENT_GRAPH

  // Find out how much free memory is on the GPU
  size_t freeMem, totalMem;
  ERR_CHK(cudaMemGetInfo(&freeMem, &totalMem));
  // Allocate 90% of it
  size_t allocMem = (freeMem * 0.9);
  ERR_CHK(cudaMalloc(&d_confVertices, allocMem));

  ERR_CHK(cudaDeviceSynchronize());

  buildCooConfGraphDevice(d_pauliEnc, pauliEncSize, d_colList, n, T, d_confOffsets, d_confVertices, d_nConflicts);
  ERR_CHK(cudaDeviceSynchronize());
  // Read d_nConflicts from GPU
  ERR_CHK(cudaMemcpy(&nConflicts, d_nConflicts, sizeof(OffsetTy), cudaMemcpyDeviceToHost));
  h_confVertices.resize(nConflicts*2);
  cubInclusiveSum((void *)d_confOffsetsCnt, n, d_confOffsets);
  cudaDeviceSynchronize();
  // if nConflicts * sizeof(NODE_T) < half of allocMem, then we can use the memory
  // allocated for d_confVertices
  if(nConflicts * 2 < (allocMem/sizeof(NODE_T))/2){
    std::cout << "Fits: " << nConflicts * 2 * sizeof(NODE_T) << " < " << allocMem/2 << std::endl;
    NODE_T *d_confCsr = d_confVertices + nConflicts*2;
    buildCsrConfGraphDevice(n, d_confOffsets, d_confOffsetsCnt, d_confVertices, d_confCsr, nConflicts);
    ERR_CHK(cudaDeviceSynchronize());
    // Read d_confCsr from GPU
    ERR_CHK(cudaMemcpy(h_confOffsets.data(), d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy), cudaMemcpyDeviceToHost));
    ERR_CHK(cudaMemcpy(h_confVertices.data(), d_confCsr, h_confVertices.size() * sizeof(NODE_T), cudaMemcpyDeviceToHost));
    ERR_CHK(cudaDeviceSynchronize());
    std::cout << h_confVertices[0] << " " << h_confVertices[nConflicts*2-1] << std::endl;
  }
  else{
    // Too Large for GPU memory, use CPU to post-process
    ERR_CHK(cudaMemcpy(h_confOffsets.data(), d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy), cudaMemcpyDeviceToHost));
    std::vector<Edge> h_cooVerticesTmp(nConflicts);
    ERR_CHK(cudaMemcpy(h_cooVerticesTmp.data(), d_confVertices, h_cooVerticesTmp.size() * sizeof(NODE_T), cudaMemcpyDeviceToHost));
    std::vector<std::atomic<OffsetTy>> h_confOffsetsTmp(n);
    for(NODE_T i = 0; i < n; i++){
      std::atomic_init(&h_confOffsetsTmp[i], 0);
    }
    #pragma omp parallel for
    for(NODE_T i = 0; i < nConflicts; i++){
      Edge e = h_cooVerticesTmp[i];
      OffsetTy offset = h_confOffsets[e.u] + std::atomic_fetch_add(&h_confOffsetsTmp[e.u], 1);
      h_confVertices[offset] = e.v;
      offset = h_confOffsets[e.v] + std::atomic_fetch_add(&h_confOffsetsTmp[e.v], 1);
      h_confVertices[offset] = e.u;
    }
  }
  std::cout << "nConflicts: " << nConflicts << std::endl;
  #pragma omp parallel for
  for(NODE_T i = 0; i < n; i++) {
    // Sort each vertex's edgelist
    std::sort(h_confVertices.begin() + h_confOffsets[i], h_confVertices.begin() + h_confOffsets[i+1]);
  }
  #else // COMPLEMENT_GRAPH
  // Compute what the offsets would be for a complement graph
  buildCooCompGraphDevice(d_pauliEnc, pauliEncSize, d_colList, n, T, d_confOffsets, d_nConflicts);
  ERR_CHK(cudaDeviceSynchronize());
  ERR_CHK(cudaMemcpy(&nConflicts, d_nConflicts, sizeof(OffsetTy), cudaMemcpyDeviceToHost));
  ERR_CHK(cudaMemcpy(h_confOffsets.data(), d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy), cudaMemcpyDeviceToHost));
  ERR_CHK(cudaDeviceSynchronize());

  // Find lowest degree and highest degree vertex
  NODE_T minDegree = n;
  NODE_T maxDegree = 0;
  double avgDegree = 0;
  for(NODE_T i = 0; i < n; i++) {
    NODE_T degree = h_confOffsets[i];
    if(degree < minDegree) {
      minDegree = degree;
    }
    if(degree > maxDegree) {
      maxDegree = degree;
    }
    avgDegree += degree;
  }
  avgDegree /= n;
  // Calculate variance
  double variance = 0;
  for(NODE_T i = 0; i < n; i++) {
    NODE_T degree = h_confOffsets[i];
    variance += (degree - avgDegree) * (degree - avgDegree);
  }
  variance /= n;
  // Print values
  std::cout << "vertices: " << n << std::endl;
  std::cout << "edges: " << nConflicts << std::endl;
  std::cout << "pauli mtx per string: " << jsongraph.getPauliLength() << std::endl;
  std::cout << "minDegree: " << minDegree << std::endl;
  std::cout << "maxDegree: " << maxDegree << std::endl;
  std::cout << "avgDegree: " << avgDegree << std::endl;
  std::cout << "variance: " << variance << std::endl;
  // Exit
  exit(0);
  #endif // COMPLEMENT_GRAPH
}

void buildConfGraphGpu (ClqPart::JsonGraph &jsongraph) {
  pauliEncSize = jsongraph.getPauliEncSize();
  std::vector<uint32_t> &h_pauliEnc = jsongraph.getEncodedData();
  ERR_CHK(cudaMalloc(&d_pauliEnc, h_pauliEnc.size() * sizeof(NODE_T)));
  ERR_CHK(cudaMemcpy(d_pauliEnc, h_pauliEnc.data(), h_pauliEnc.size() * sizeof(NODE_T), cudaMemcpyHostToDevice));

  h_confOffsets.resize(n);
  ERR_CHK(cudaMalloc(&d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy)));
  #ifndef COMPLEMENT_GRAPH
  h_confVertices.resize(n * n);
  ERR_CHK(cudaMalloc(&d_confVertices, h_confVertices.size() * sizeof(NODE_T)));
  #endif // !COMPLEMENT_GRAPH

  ERR_CHK(cudaDeviceSynchronize());

  // Create d_nConflicts and initialize to 0
  OffsetTy *d_nConflicts;
  ERR_CHK(cudaMalloc(&d_nConflicts, sizeof(OffsetTy)));
  ERR_CHK(cudaMemset(d_nConflicts, 0, sizeof(OffsetTy)));

  // Call function to send to GPU to build conf graph
  #ifdef COMPLEMENT_GRAPH
  buildCompGraphDevice(d_pauliEnc, pauliEncSize, d_colList, n, T, d_confOffsets, d_confVertices, d_nConflicts);
  #else // !COMPLEMENT_GRAPH
  buildConfGraphDevice(d_pauliEnc, pauliEncSize, d_colList, n, T, d_confOffsets, d_confVertices, d_nConflicts);
  #endif // !COMPLEMENT_GRAPH
  ERR_CHK(cudaDeviceSynchronize());
  // Read d_nConflicts from GPU
  ERR_CHK(cudaMemcpy(&nConflicts, d_nConflicts, sizeof(OffsetTy), cudaMemcpyDeviceToHost));
  // Print nConflicts
  nConflicts /= 2;
  // Copy h_confOffsets and h_confVertices from GPU
  ERR_CHK(cudaMemcpy(h_confOffsets.data(), d_confOffsets, h_confOffsets.size() * sizeof(OffsetTy), cudaMemcpyDeviceToHost));
  ERR_CHK(cudaMemcpy(h_confVertices.data(), d_confVertices, h_confVertices.size() * sizeof(NODE_T), cudaMemcpyDeviceToHost));

  ERR_CHK(cudaDeviceSynchronize());
  std::cout << "nConflicts: " << nConflicts << std::endl;
  // Exit program
  for(NODE_T i = 0; i < n; i++) {
    confAdjList[i].reserve(h_confOffsets[i]);
    for(NODE_T j = 0; j < h_confOffsets[i]; j++) {
      confAdjList[i].push_back(h_confVertices[i * n + j]);
    }
    std::sort(confAdjList[i].begin(), confAdjList[i].end());
  }
  #ifdef COMPLEMENT_GRAPH
  // Find lowest degree and highest degree vertex
  NODE_T minDegree = n;
  NODE_T maxDegree = 0;
  double avgDegree = 0;
  for(NODE_T i = 0; i < n; i++) {
    if(h_confOffsets[i] < minDegree) {
      minDegree = h_confOffsets[i];
    }
    if(h_confOffsets[i] > maxDegree) {
      maxDegree = h_confOffsets[i];
    }
    avgDegree += h_confOffsets[i];
  }
  avgDegree /= n;
  // Calculate variance
  double variance = 0;
  for(NODE_T i = 0; i < n; i++) {
    variance += (h_confOffsets[i] - avgDegree) * (h_confOffsets[i] - avgDegree);
  }
  variance /= n;
  // Print values
  std::cout << "vertices: " << n << std::endl;
  std::cout << "edges: " << nConflicts << std::endl;
  std::cout << "pauli mtx per string: " << jsongraph.getPauliLength() << std::endl;
  std::cout << "minDegree: " << minDegree << std::endl;
  std::cout << "maxDegree: " << maxDegree << std::endl;
  std::cout << "avgDegree: " << avgDegree << std::endl;
  std::cout << "variance: " << variance << std::endl;
  // Exit
  exit(0);
  #endif // COMPLEMENT_GRAPH
}

void confColorGreedyCSR() {
  
  double t1 = omp_get_wtime();
  std::cout<<"# of conflicting edges: "<<nConflicts<<std::endl; 
  //Buckets of size T;
  std::vector<std::vector<NODE_T> > verBucket(T+1);
  //to stoe the position of vertex in the corresponding bucket
  std::vector<NODE_T> verLocation(n);
  
  NODE_T vMin = T; 

  std::mt19937 engine(2034587);
  std::uniform_int_distribution<NODE_T> uniform_dist(0, colThreshold-1);

  //# of verties processed. For stopping condition later
  NODE_T vtxProcessed = 0;
  for(NODE_T i=0;i<n;i++) {
    //if this is a non-conflicting vertex we can color it arbitrarily from the 
    //list of colors.
    if(h_confOffsets[i+1] - h_confOffsets[i] == 0) {
      //std::cout<<"coloring non-conflicting edge"<<std::endl;
      uniform_dist = std::uniform_int_distribution<NODE_T>(0, 
          colList[i].size()-1);
      colors[i] = colList[i][uniform_dist(engine)];
      vtxProcessed++;
      continue;
    }
    //otherwise place it in appropriate bucket
    NODE_T t = colList[i].size();
    verBucket[t-1].push_back(i); 
    verLocation[i] = verBucket[t-1].size()-1;
    
    //the minimum bucket length
    if(t < vMin) {
      vMin = t; 
    }
  }
 
  //keep processing until we see all the vertices. 
  while (vtxProcessed < n) {
    for(NODE_T i = vMin-1; i < T;i++) {
      //if this bucket has elements
      //select a random vertex and swap it with the last element
      //to efficiently remove this vertex from the bucket.
      //std::cout<<"Bucket: "<<i<<"size: "<<verBucket[i].size()<<std::endl;
      if(verBucket[i].empty() == false) {
        uniform_dist = std::uniform_int_distribution<NODE_T>(0, 
            verBucket[i].size()-1);
        NODE_T selectedVtxLoc = uniform_dist(engine);
        NODE_T selectedVtx = verBucket[i].at(selectedVtxLoc);

        //update the verLocation array of the last element
        verLocation[verBucket[i].back()] = selectedVtxLoc;
        
        //swap the vertex with the last element
        std::swap(verBucket[i][selectedVtxLoc],verBucket[i].back());
        //remove the vertex
        verBucket[i].pop_back();

        //attempt to color this vertex
        NODE_T colInd = attemptToColor(selectedVtx);
        if(colInd>=0) {
          //remove col at colInd of the vtx is not necessary
          //std::swap(colList[selectedVtx][colInd],colList[selectedVtx].back());
          //colList[selectedVtx].pop_back();

          fixBucketsCSR(selectedVtx,colors[selectedVtx],vtxProcessed,verBucket,verLocation,vMin); 
          //std::cout<<"updated VMin: "<<vMin<<std::endl;
        }
        else {
          colors[selectedVtx] = -2;
          invalidVertices.push_back(selectedVtx); 
        }
        vtxProcessed++;
        break;
      } 
    }
  }
  std::cout<<"# of invalid vertices: "<<invalidVertices.size()
    <<std::endl;

 /* 
  NODE_T invalid = 0;
  for(NODE_T i=0;i<n;i++) {
    if(colors[i] == -2)
     invalid++; 
  }
  std::cout<<invalid<<std::endl;
  */

  confColorTime = omp_get_wtime() - t1;
  
  std::cout<<"Conflict Coloring Time: "<<confColorTime<<std::endl;

  nColors = *std::max_element(colors.begin(),colors.end()) + 1;

}

void fixBucketsCSR(NODE_T vtx, NODE_T col, NODE_T &vtxProcessed, 
    std::vector< std::vector<NODE_T> > &verBucket, std::vector<NODE_T> &verLoc, NODE_T &vMin) {
  for(auto v_id = h_confOffsets[vtx]; v_id < h_confOffsets[vtx+1]; v_id++) {
    auto v = h_confVertices[v_id];
    if(colors[v] == -1) {
      auto it = std::find(colList[v].begin(),colList[v].end(),col);
      if(it != colList[v].end()) {
        //fix the location of the last vertex 
        verLoc[verBucket[colList[v].size()-1].back()] = verLoc[v];
        //now swap v with last element
        std::swap(verBucket[colList[v].size()-1][verLoc[v]],verBucket[colList[v].size()-1].back());
        //remove v from this bucket
        verBucket[colList[v].size()-1].pop_back();
        
        //swap the colr with last color
        std::swap(*it,colList[v].back());
        //remove the color
        colList[v].pop_back(); 
        if(colList[v].empty() == true) {
          vtxProcessed++;
          invalidVertices.push_back(v);
          colors[v] = -2;
        }
        else {
          if (colList[v].size() < vMin) 
            vMin = colList[v].size();

          verBucket[colList[v].size()-1].push_back(v); 
          verLoc[v] = verBucket[colList[v].size()-1].size()-1;
        }
      }
      
    }
  } 

  
}
#endif // ENABLE_GPU

//This function sort the vertices of the conflict graph w.r.t to their degrees.
//The idea is to color the vertex with highest degree first, since this represents
//the most conflicted with other vertices. Does not perform that well in practice.
//We are not using it at the moment.
void orderConfVertices() {
   
  std::vector<NODE_T > degreeConf(n,0);

  for(NODE_T i=0;i<n;i++) {
    degreeConf[i]= confAdjList[i].size(); 
  }
  std::iota(vertexOrder.begin(),vertexOrder.end(),0);

  std::stable_sort(vertexOrder.begin(),vertexOrder.end(),
      [&degreeConf] (NODE_T t1, NODE_T t2) {return degreeConf[t1] > degreeConf[t2];});


}

//This function attempt to color the conflicting graph with largest degree heuristics.
//Does not perform well
void confColor() {
  
  std::cout<<"conflicting edges: "<<nConflicts<<std::endl; 
  orderConfVertices();
  for(NODE_T i:vertexOrder) {
    if(confAdjList[i].empty() == false) {
      for(auto col:colList[i]) {
        bool flag = true;
        for(auto v:confAdjList[i]) {
          if (colors[v] == col) {
            flag = false;
            break; 
          }
        }     

        if(flag == true) {
          colors[i] = col;
          break; 
        }
      }  
      if(colors[i] == -1)
        invalidVertices.push_back(i);
    } 
    else
      colors[i] = colList[i][0];
    //std::cout<<i<<" "<<colors[i]<<std::endl;
  }
  std::cout<<"# of vertices can not be colored: "<<invalidVertices.size()
    <<std::endl;
  nColors = *std::max_element(colors.begin(),colors.end()) + 1;
}

void confColorGreedy() {
  
  double t1 = omp_get_wtime();
  std::cout<<"# of conflicting edges: "<<nConflicts<<std::endl; 
  //Buckets of size T;
  std::vector<std::vector<NODE_T> > verBucket(T+1);
  //to stoe the position of vertex in the corresponding bucket
  std::vector<NODE_T> verLocation(n);
  
  NODE_T vMin = T; 

  std::mt19937 engine(2034587);
  std::uniform_int_distribution<NODE_T> uniform_dist(0, colThreshold-1);

  //# of verties processed. For stopping condition later
  NODE_T vtxProcessed = 0;
  for(NODE_T i=0;i<n;i++) {
    //if this is a non-conflicting vertex we can color it arbitrarily from the 
    //list of colors.
    if(confAdjList[i].empty() == true) {
      //std::cout<<"coloring non-conflicting edge"<<std::endl;
      uniform_dist = std::uniform_int_distribution<NODE_T>(0, 
          colList[i].size()-1);
      colors[i] = colList[i][uniform_dist(engine)];
      vtxProcessed++;
      continue;
    }
    //otherwise place it in appropriate bucket
    NODE_T t = colList[i].size();
    verBucket[t-1].push_back(i); 
    verLocation[i] = verBucket[t-1].size()-1;
    
    //the minimum bucket length
    if(t < vMin) {
      vMin = t; 
    }
  }
 
  //keep processing until we see all the vertices. 
  while (vtxProcessed < n) {
    for(NODE_T i = vMin-1; i < T;i++) {
      //if this bucket has elements
      //select a random vertex and swap it with the last element
      //to efficiently remove this vertex from the bucket.
      //std::cout<<"Bucket: "<<i<<"size: "<<verBucket[i].size()<<std::endl;
      if(verBucket[i].empty() == false) {
        uniform_dist = std::uniform_int_distribution<NODE_T>(0, 
            verBucket[i].size()-1);
        NODE_T selectedVtxLoc = uniform_dist(engine);
        NODE_T selectedVtx = verBucket[i].at(selectedVtxLoc);

        //update the verLocation array of the last element
        verLocation[verBucket[i].back()] = selectedVtxLoc;
        
        //swap the vertex with the last element
        std::swap(verBucket[i][selectedVtxLoc],verBucket[i].back());
        //remove the vertex
        verBucket[i].pop_back();

        //attempt to color this vertex
        NODE_T colInd = attemptToColor(selectedVtx);
        if(colInd>=0) {
          //remove col at colInd of the vtx is not necessary
          //std::swap(colList[selectedVtx][colInd],colList[selectedVtx].back());
          //colList[selectedVtx].pop_back();

          fixBuckets(selectedVtx,colors[selectedVtx],vtxProcessed,verBucket,verLocation,vMin); 
          //std::cout<<"updated VMin: "<<vMin<<std::endl;
        }
        else {
          colors[selectedVtx] = -2;
          invalidVertices.push_back(selectedVtx); 
        }
        vtxProcessed++;
        break;
      } 
    }
  }
  std::cout<<"# of invalid vertices: "<<invalidVertices.size()
    <<std::endl;

 /* 
  NODE_T invalid = 0;
  for(NODE_T i=0;i<n;i++) {
    if(colors[i] == -2)
     invalid++; 
  }
  std::cout<<invalid<<std::endl;
  */

  confColorTime = omp_get_wtime() - t1;
  
  std::cout<<"Conflict Coloring Time: "<<confColorTime<<std::endl;

  nColors = *std::max_element(colors.begin(),colors.end()) + 1;

}

// #endif

  // #ifdef ENABLE_GPU
  // void buildConfGraphGpu(ClqPart::JsonGraph &jsongraph);
  // #endif
  // void buildConfGraph( ClqPart::JsonGraph &);
  // void confColor();
  // void confColorGreedy();
  // void naiveGreedyColor(std::vector<NODE_T> vertList, ClqPart::JsonGraph &jsongraph,NODE_T offset);
  // void confColorRand();
  // void orderConfVertices();

  std::vector< std::vector<NODE_T> >& getConfAdjList() { return confAdjList; }
  std::vector<NODE_T>& getConfVertices() { return confVertices; }
  std::vector<NODE_T>& getInvVertices() { return invalidVertices; }
  std::vector<NODE_T> getColors() { return colors; }
  NODE_T getNumColors() {return nColors;}

private:
  //This function assign random list of colors from the Palette.
void assignListColor() {

  std::mt19937 engine(2034587);
  std::uniform_int_distribution<NODE_T> uniform_dist(0, colThreshold-1);
  
  //We want the colors to be not repeating. Thus initially the list size is
  //same for each vertex
  std::vector<bool> isPresent(colThreshold);

  std::cout<<"Assigning random list color " << "Palette Size: "<<colThreshold
            <<" "<<"List size: "<<T<<std::endl;
  
  double t1 = omp_get_wtime();
  #ifdef ENABLE_GPU
  h_colList.reserve(n * T);
  #endif // ENABLE_GPU
  for (NODE_T i= 0; i<n ; i++) {
    std::fill(isPresent.begin(),isPresent.end(),false);
    for (NODE_T j=0; j<T ; j++) {
      NODE_T col;
      do {
        col = uniform_dist(engine);
      }while(isPresent[col] == true);

      colList[i].push_back(col); 
      isPresent[col] = true;
    } 
    std::sort(colList[i].begin(),colList[i].end());
    #ifdef ENABLE_GPU
    h_colList.insert(h_colList.end(), colList[i].begin(), colList[i].end());
    #endif
  }
  #ifdef ENABLE_GPU
  // Copy h_colList to GPU
  cudaError_t err;
  ERR_CHK(cudaMalloc(&d_colList, h_colList.size() * sizeof(NODE_T)));
  ERR_CHK(cudaMemcpy(d_colList, h_colList.data(), h_colList.size() * sizeof(NODE_T), cudaMemcpyHostToDevice));
  #endif // ENABLE_GPU
  assignTime = omp_get_wtime() - t1;
  std::cout<<"Assignment Time: "<<assignTime<<std::endl;

}
  /********************************************************************************
//The next three functions are not necessary at this moment. They are used to color
//the invalid vertices. But since we do not know the subgraph induced by the invalid
//vertices coloring them only with conflict graph is wrong.  
void populateCandColors( NODE_T u, std::vector<NODE_T> &candColors) {
  for(auto v:confAdjList[u] ) {
    if (confColors[v] !=  -1) {
      candColors[ confColors[v] ] = u; 
    } 
  }
}

NODE_T firstAvailColor( NODE_T u, std::vector<NODE_T> &candColors) {
  for(NODE_T v=0;v < candColors.size(); v++ ) {
    if ( candColors[v] != u) {
      return v; 
    } 
  }
  return -1;
}
void greedyColor(NODE_T offset) {
  std::vector<NODE_T> candColors(n,-1);
 
  for(auto u:invalidVertices ) {
    //std::cout<<u<<std::endl;
    populateCandColors(u,candColors); 
    
    confColors[u] = firstAvailColor(u,candColors);
    colors[u] = confColors[u]+offset;
  }
}
********************************************************************************/
  //The function given a vtx, properly color the vertx using a random color from the list.
//the coloring is guaranteed but I first name it attempt. kept in this way.
NODE_T attemptToColor(NODE_T vtx) {

  std::mt19937 engine(213857);
  std::uniform_int_distribution<NODE_T> uniform_dist(0,colList[vtx].size()-1);

  auto colInd = uniform_dist(engine);
  
  auto col = colList[vtx].at(colInd);
  colors[vtx] = col;
  return colInd;
}
  //The next two functions fixBuckets and ConfColorGreedy are needed for coloring
//the conflict graph. This heuristic color the vertices with the smallest list size
//first, since they have less options. Works well in practice. 
void fixBuckets(NODE_T vtx, NODE_T col, NODE_T &vtxProcessed, 
    std::vector< std::vector<NODE_T> > &verBucket, std::vector<NODE_T> &verLoc, NODE_T &vMin) {
  for( auto v:confAdjList[vtx] ) {
    if(colors[v] == -1) {
      auto it = std::find(colList[v].begin(),colList[v].end(),col);
      if(it != colList[v].end()) {
        //fix the location of the last vertex 
        verLoc[verBucket[colList[v].size()-1].back()] = verLoc[v];
        //now swap v with last element
        std::swap(verBucket[colList[v].size()-1][verLoc[v]],verBucket[colList[v].size()-1].back());
        //remove v from this bucket
        verBucket[colList[v].size()-1].pop_back();
        
        //swap the colr with last color
        std::swap(*it,colList[v].back());
        //remove the color
        colList[v].pop_back(); 
        if(colList[v].empty() == true) {
          vtxProcessed++;
          invalidVertices.push_back(v);
          colors[v] = -2;
        }
        else {
          if (colList[v].size() < vMin) 
            vMin = colList[v].size();

          verBucket[colList[v].size()-1].push_back(v); 
          verLoc[v] = verBucket[colList[v].size()-1].size()-1;
        }
      }
      
    }
  } 

  
}

};
