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

#ifdef ENABLE_GPU
// #include <cuda.h>
#include <cuda_runtime.h>
#endif

#include "ClqPart/graph.h" 
#include "ClqPart/JsonGraph.h"
#include "ClqPart/cuPaletteCol.cuh"
#include <random>

#include <omp.h>

bool findFirstCommonElement(const std::vector<NODE_T>& vec1, const std::vector<NODE_T>& vec2);

class PaletteColor {
  
  NODE_T n;
  NODE_T colThreshold;
  std::vector<std::vector<NODE_T> > colList;
  std::vector<NODE_T> colors;
  std::vector<NODE_T> confColors;
  NODE_T nColors;
  NODE_T nConflicts;

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
  std::vector<NODE_T> h_confOffsets;
  NODE_T *d_confOffsets;
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

  void buildStreamConfGraph( NODE_T u, NODE_T v ); 
  #ifdef ENABLE_GPU
  void buildConfGraphGpu(ClqPart::JsonGraph &jsongraph);
  #endif
  // void buildConfGraph( ClqPart::JsonGraph &);
  void confColor();
  void confColorGreedy();
  // void naiveGreedyColor(std::vector<NODE_T> vertList, ClqPart::JsonGraph &jsongraph,NODE_T offset);
  void confColorRand();
  void orderConfVertices();

  std::vector< std::vector<NODE_T> >& getConfAdjList() { return confAdjList; }
  std::vector<NODE_T>& getConfVertices() { return confVertices; }
  std::vector<NODE_T>& getInvVertices() { return invalidVertices; }
  std::vector<NODE_T> getColors() { return colors; }
  NODE_T getNumColors() {return nColors;}

private:
  void assignListColor();
  void populateCandColors(NODE_T, std::vector<NODE_T> &);
  void greedyColor(NODE_T);
  NODE_T firstAvailColor(NODE_T, std::vector<NODE_T> &);
  NODE_T attemptToColor(NODE_T);
  void fixBuckets(NODE_T , NODE_T , NODE_T &, std::vector< std::vector<NODE_T> > &
      , std::vector<NODE_T> &, NODE_T &); 

};
