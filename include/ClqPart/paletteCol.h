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
#include "ClqPart/JsonGraph.h"
#include <random>

#include <omp.h>

struct PalColStat {
  NODE_T n;
  EDGE_T m;
  EDGE_T mConf;
  NODE_T palSz;
  NODE_T lstSz;
  NODE_T nColors;
  double assignTime;
  double confBuildTime;
  double confColorTime; 
  double invColorTime;
};

class PaletteColor {
  
  NODE_T n;
  NODE_T colThreshold;
  NODE_T T;
  std::vector<std::vector<NODE_T> > colList;
  std::vector<NODE_T> colors;
  NODE_T nColors;
  EDGE_T nConflicts;
  int level;
  
  std::vector<NODE_T> invalidVertices;
  std::vector<NODE_T> vertexOrder;
  std::vector<std::vector< NODE_T> > confAdjList;
  std::vector<PalColStat> palStat;

public:
  PaletteColor( NODE_T n1, NODE_T colThresh, float alpha=1, NODE_T lst_sz = -1) {
    n = n1;
    colThreshold = colThresh; 
    colors.resize(n,-1);
    nColors = 0;
    nConflicts = 0;
    if(lst_sz < 0)
      T =  static_cast<NODE_T> (alpha*log(n));
    else
      T = lst_sz;

    if(T>colThresh)
      T = colThresh;

    confAdjList.resize(n);
    colList.resize(n);
    invalidVertices.clear();
    palStat.push_back({n,-1,-1,colThreshold,T,0,0.0,0.0,0.0,0.0});
    level = 0;
    assignListColor();
  }
  //initialize for the recursive implementation. 
  void reInit(std::vector<NODE_T> & nodeList,NODE_T colThresh, float alpha=1, NODE_T lst_sz = -1) {
    colThreshold = colThresh; 
    nConflicts = 0;
    //The nodes in the node List need to have color -1.
    for(auto u:nodeList) {
      colors[u] = -1; 
      colList[u].clear();
      confAdjList[u].clear();
    }
    if(lst_sz < 0)
      T =  static_cast<NODE_T> (alpha*log(nodeList.size()));
    else
      T = lst_sz;

    if(T>colThresh)
      T = colThresh;
    invalidVertices.clear();

    palStat.push_back({nodeList.size(),-1,-1,colThreshold,T,0,0.0,0.0,0.0});
    level = level + 1;
    assignListColor(nodeList,getNumColors());
  }

  void buildStreamConfGraph( NODE_T u, NODE_T v ); 
  void buildConfGraph( ClqPart::JsonGraph &);
  void buildConfGraph( ClqPart::JsonGraph &,std::vector<NODE_T> &);

  void confColor();
  void confColorGreedy(std::vector<NODE_T> &);
  void confColorGreedy();

  void naiveGreedyColor(std::vector<NODE_T> vertList, ClqPart::JsonGraph &jsongraph,NODE_T offset);
  void confColorRand();
  void orderConfVertices();
  bool checkValidity(ClqPart::JsonGraph &);

  std::vector< std::vector<NODE_T> >& getConfAdjList() { return confAdjList; }
  //std::vector<NODE_T>& getConfVertices() { return confVertices; }
  std::vector<NODE_T>& getInvVertices() { return invalidVertices; }
  std::vector<NODE_T> getColors() { return colors; }
  NODE_T getNumColors() {return nColors;}
  EDGE_T getNumConflicts() {return nConflicts;}
  PalColStat getPalStat(int lev=0) {return palStat[lev];}

private:
  void assignListColor();
  void assignListColor(std::vector<NODE_T> &, NODE_T);
  void populateCandColors(NODE_T, std::vector<NODE_T> &);
  void greedyColor(NODE_T);
  NODE_T firstAvailColor(NODE_T, std::vector<NODE_T> &);
  NODE_T attemptToColor(NODE_T);
  void fixBuckets(NODE_T , NODE_T , NODE_T &, std::vector< std::vector<NODE_T> > &
      , std::vector<NODE_T> &, NODE_T &); 

};
