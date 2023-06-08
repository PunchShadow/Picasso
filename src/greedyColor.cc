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

#include <numeric>
#include <iostream>

#include "ClqPart/greedyColor.h"
#include "ClqPart/utility.h"



void GreedyColor::populateCandColors( NODE_T u, std::vector<NODE_T> &candColors) {
  for( auto i=G.IA[u] ; i<G.IA[u+1] ; i++ ) {
    if (colors[G.JA[i]] !=  -1) {
      candColors[ colors[G.JA[i]] ] = u; 
    } 
  }
}

NODE_T GreedyColor::firstAvailColor( NODE_T u, std::vector<NODE_T> &candColors) {
  for(NODE_T v=0;v < candColors.size(); v++ ) {
    if ( candColors[v] != u) {
      return v; 
    } 
  }
  return -1;
}

void GreedyColor::orderLargestFirst() {
  NODE_T n = G.numberOfNodes();
  
  std::vector<NODE_T> degrees(n,0);

  getDegrees(G,degrees);
 
  if(vertexOrder.size() != degrees.size()) { 
    vertexOrder.resize(degrees.size());
  }
  std::iota(vertexOrder.begin(),vertexOrder.end(),0);

  std::stable_sort(vertexOrder.begin(),vertexOrder.end(),
      [&degrees] (NODE_T t1, NODE_T t2) {return degrees[t1] > degrees[t2];});

  /*for(auto v:vertexOrder ) {
    std::cout<<v<<"\n"; 
  }*/
}

void GreedyColor::colorVertices() {

  double t1 = omp_get_wtime();
  NODE_T n = G.numberOfNodes();
  std::vector<NODE_T> candColors(n,-1);
 
  for(auto u:vertexOrder ) {
      populateCandColors(u,candColors); 
      colors[u] = firstAvailColor(u,candColors);
  }
  nColors = *std::max_element(colors.begin(),colors.end()) + 1;
  colorTime = omp_get_wtime() - t1;
}

void GreedyColor::orderVertices() {

  double t1 = omp_get_wtime();
  if (orderName == "LARGEST_FIRST" ) {
    orderLargestFirst(); 
  }
  orderTime = omp_get_wtime() - t1;
}
