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

#include "ClqPart/graph.h" 
#include <random>

#include <omp.h>

class PaletteColor {
  
  NODE_T n;
  NODE_T colThreshold;
  std::vector<std::vector<NODE_T> > colList;
  std::vector<NODE_T> colors;
  std::vector<NODE_T> confColors;
  NODE_T nColors;
  NODE_T nConflicts;
  
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
  void buildStreamConfGraph( NODE_T u, NODE_T v ); 
  void confColor();
  void confColorGreedy();
  void confColorRand();
  void orderConfVertices();

  std::vector< std::vector<NODE_T> >& getConfAdjList() { return confAdjList; }
  std::vector<NODE_T>& getConfVertices() { return confVertices; }
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
