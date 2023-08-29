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
  
  cxxopts::Options options("palettecol", "read json pauli string files and color the graph using palette coloring algorithm"); 
  options.add_options()
    ("in,infile", "input json file name", cxxopts::value<std::string>())
    ("t,target", "target color", cxxopts::value<NODE_T>())
    ("a,alpha", "coefficient to log(n)", cxxopts::value<float>()->default_value("1.0"))
    ("s,stream", "streaming construction", cxxopts::value<bool>()->default_value("false"))
    ("l,list", "list size ", cxxopts::value<NODE_T>()->default_value("-1"))
    ("h,help", "print usage")
    ;

  std::string inFname;
  NODE_T target,list_size;
  float alpha;
  bool isStream;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    target = result["target"].as<NODE_T>();
    alpha = result["alpha"].as<float>();
    isStream = result["stream"].as<bool>();
    list_size = result["list"].as<NODE_T>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  ClqPart::JsonGraph jsongraph(inFname,isStream); 
  NODE_T n = jsongraph.numOfData();
  if(list_size >=0)
    std::cout<<"Since list size is given, ignoring alpha"<<std::endl;
  PaletteColor palcol(n,target,alpha,list_size);

  double t1 = omp_get_wtime();
  if(isStream){
    ClqPart::Edge e;
    while(jsongraph.nextEdge(e)) {
       palcol.buildStreamConfGraph(e.u,e.v); 
    
    }
  }
  else { 
    palcol.buildConfGraph(jsongraph);
  }
  double createConfTime = omp_get_wtime() - t1;
  std::cout<<"Conflict graph construction time: "<<createConfTime<<std::endl;

  std::vector< std::vector<NODE_T> > confEdges = palcol.getConfAdjList();
  palcol.confColorGreedy();
  std::vector<NODE_T> colors = palcol.getColors();
  std::cout<<"# of colors: " <<palcol.getNumColors()<<std::endl;
  
  std::vector <NODE_T>  invVert = palcol.getInvVertices();
  if(invVert.empty() == false) {
    palcol.naiveGreedyColor(invVert, jsongraph, palcol.getNumColors());
    std::cout<<"# of Final colors: " <<palcol.getNumColors()<<std::endl;
  }
  return 0;
}  
