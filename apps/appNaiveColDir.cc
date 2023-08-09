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

#include <iostream>
#include <cstring>
#include<fstream>
#include <random>

#include "ClqPart/JsonGraph.h"
#include "ClqPart/paletteCol.h"
#include "ClqPart/input.h"
#include "cxxopts/cxxopts.hpp"
#include "ClqPart/utility.h"
#include "ClqPart/greedyColor.h"


int main(int argC, char *argV[]) {
  
  cxxopts::Options options("naivecol", "read json pauli string files and color the graph using naive coloring algorithm"); 
  options.add_options()
    ("in,infile", "input json file name", cxxopts::value<std::string>())
//("out,outfile", "output color file name", cxxopts::value<std::string>()->default_value(""))
    //("o,order", "LARGEST_FIRST,SMALLEST_LAST,NATURAL,RANDOM,DYNAMIC_LARGEST_FIRST,INCIDENCE_DEGREE",cxxopts::value<std::string>()->default_value("LARGEST_FIRST"))
    ("h,help", "print usage")
    ;

  //std::string inFname, outFname,order;
  std::string inFname;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  ClqPart::JsonGraph jsongraph(inFname,true); 
  NODE_T n = jsongraph.numOfData();
  //jsongraph.printData();

  std::vector<NODE_T> forbiddenCol(n,-1);
  std::vector<NODE_T> colors(n,-1);
  
  colors[0] = 0;
  for(auto u=1 ; u<n; u++) {
    //populate the forbidden colors for u
    for(auto v=0; v<u; v++) {
      if(jsongraph.is_an_edge(v,u) == false) { 
        if (colors[v] != -1) {
          forbiddenCol[colors[v]] = u;
        }
      }
    } 
    //color u with first available color
    for( auto i=0;i<n;i++) {
      if(forbiddenCol[i] == -1) {
        colors[u] = i; 
        break;
      } 
    }
  }
  
  NODE_T nColors = *std::max_element(colors.begin(),colors.end()) + 1;
  std::cout<<jsongraph.numOfData()<<" "
    <<jsongraph.getNumEdge()<<" "
    <<nColors<<std::endl;
  /*PaletteColor palcol(n,target,alpha);

  ClqPart::Edge e;
  double t1 = omp_get_wtime();
  while(jsongraph.nextEdge(e)) {
     palcol.buildStreamConfGraph(e.u,e.v); 
  
  }
  double createConfTime = omp_get_wtime() - t1;
  std::cout<<"creating conflict graph time: "<<createConfTime<<std::endl;

  std::vector< std::vector<NODE_T> > confEdges = palcol.getConfAdjList();
  palcol.confColorGreedy();
  std::vector<NODE_T> colors = palcol.getColors();
  */
  /*
  std::vector<NODE_T> colHist(n/7,0);
  for(auto u:colors) {
    if(u>=0) colHist[u]++; 
  }
  for( auto c:colHist) {
    std::cout<<c<<std::endl;
    if(c==0)
     unassigned++; 
  }
  */
  return 0;
}  
