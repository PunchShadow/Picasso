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

void printStat( int level, PalColStat &palStat) {
  
  std::cout<<"***********Level "<<level<<"*******"<<std::endl;
  std::cout<<"Num Nodes: "<<palStat.n<<"\n";
  std::cout<<"Num Edges: "<<palStat.m<<"\n";
  std::cout<<"Avg. Deg.: "<<(double) 2*palStat.m/palStat.n<<"\n";
  std::cout<<"Palette Size: "<<palStat.palSz<<"\n";
  std::cout<<"List Size: "<<palStat.lstSz<<"\n";
  std::cout<<"Num Conflict Edges: "<<palStat.mConf<<"\n";
  std::cout<<"Conflict to Edge (%): "<<(double)palStat.mConf/palStat.m*100<<"\n";
  std::cout<<"Num Colors: " <<palStat.nColors<<std::endl;
  std::cout<<"Assign Time: "<<palStat.assignTime<<"\n";
  std::cout<<"Conf. Build Time: "<<palStat.confBuildTime<<"\n";
  std::cout<<"Conf. Color Time: "<<palStat.confColorTime<<"\n";
  std::cout<<"\n";

}

int main(int argC, char *argV[]) {
  
  cxxopts::Options options("palettecol", "read json pauli string files and color the graph using palette coloring algorithm"); 
  options.add_options()
    ("in,infile", "json file containing the pauli strings", cxxopts::value<std::string>())
    ("t,target", "palette size", cxxopts::value<NODE_T>())
    ("a,alpha", "coefficient to log(n) for list size", cxxopts::value<float>()->default_value("1.0"))
    //("s,stream", "use streaming construction", cxxopts::value<bool>()->default_value("false"))
    ("l,list", "use explicit list size", cxxopts::value<NODE_T>()->default_value("-1"))
    ("c,check", "check validity of coloring", cxxopts::value<bool>()->default_value("false"))
    ("r,recurse", "use recursive coloring", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "print usage")
    ;

  std::string inFname;
  NODE_T target,list_size;
  float alpha;
  bool isValid,isRec;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    target = result["target"].as<NODE_T>();
    alpha = result["alpha"].as<float>();
    //isStream = result["stream"].as<bool>();
    isValid = result["check"].as<bool>();
    isRec = result["recurse"].as<bool>();
    list_size = result["list"].as<NODE_T>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  ClqPart::JsonGraph jsongraph(inFname); 
  NODE_T n = jsongraph.numOfData();
  if(list_size >=0)
    std::cout<<"Since list size is given, ignoring alpha"<<std::endl;
  PaletteColor palcol(n,target,alpha,list_size);
  
  int level = 0;
  //level 0 
  //**************Conflict graph construction*****
  //double t1 = omp_get_wtime();
  /*if(isStream){
    ClqPart::Edge e;
    while(jsongraph.nextEdge(e)) {
       palcol.buildStreamConfGraph(e.u,e.v); 
    
    }
  }*/
  palcol.buildConfGraph(jsongraph);
  
  //double createConfTime = omp_get_wtime() - t1;
  //std::cout<<"Conflict graph construction time: "<<createConfTime<<std::endl;

  //std::vector< std::vector<NODE_T> > confEdges = palcol.getConfAdjList();
  palcol.confColorGreedy();
  //std::vector<NODE_T> colors = palcol.getColors();
  
  std::vector <NODE_T>  invVert = palcol.getInvVertices();
  PalColStat palStat = palcol.getPalStat(level); 
  palStat.m = jsongraph.getNumEdge();
  //palStat.mConf = palcol.getNumConflicts();
  palStat.nColors = palcol.getNumColors(); 
  printStat(level,palStat);

  if (isRec == true) {
    while(invVert.size() > 100) { 
      jsongraph.resetNumEdge();
      level++;
      if(invVert.empty() == false) {
        if(invVert.size() > 40000) alpha = 3;
        else if (invVert.size() > 20000) alpha = 2; 
        else if (invVert.size() > 5000) alpha = 1.5; 
        else alpha = 1;
        palcol.reInit(invVert,invVert.size()/8,alpha);
        palcol.buildConfGraph(jsongraph,invVert);
        palcol.confColorGreedy(invVert);
      }
      invVert = palcol.getInvVertices();
      palStat = palcol.getPalStat(level); 
      palStat.m = jsongraph.getNumEdge();
      //palStat.mConf = palcol.getNumConflicts();
      palStat.nColors = palcol.getNumColors(); 
      printStat(level,palStat);
    }
  }
  if(invVert.empty() == false) {
    std::cout<<"Final Num invalid Vert: "<<invVert.size()<<"\n";
    palcol.naiveGreedyColor(invVert, jsongraph, palcol.getNumColors());
    palStat = palcol.getPalStat(level); 
    std::cout<<"Naive Color TIme: "<<palStat.invColorTime<<"\n";
  }
  
  std::cout<<"# of Final colors: " <<palcol.getNumColors()<<std::endl;
  
  if(isValid) {
    if(palcol.checkValidity(jsongraph)) 
      std::cout<<"Coloring valid"<<std::endl;
    else
      std::cout<<"Coloring invalid"<<std::endl;
  }

  return 0;
}  
