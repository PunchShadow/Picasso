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

#include "ClqPart/paletteCol.h"
#include "ClqPart/input.h"
#include "cxxopts/cxxopts.hpp"
#include "ClqPart/utility.h"

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
  
  cxxopts::Options options("palettecol", "read mtx and color the graph using palette coloring algorithm"); 
  options.add_options()
    ("in,infile", "input mtx file name", cxxopts::value<std::string>())
    ("t,target", "target color", cxxopts::value<NODE_T>())
    ("a,alpha", "coefficient to log(n)", cxxopts::value<float>()->default_value("1.0"))
    ("l,list", "list size ", cxxopts::value<NODE_T>()->default_value("-1"))
    ("c,check", "check validity of coloring", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "print usage")
    ;

  std::string inFname;
  NODE_T target,list_size;
  float alpha;
  bool isValid;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    target = result["target"].as<NODE_T>();
    alpha = result["alpha"].as<float>();
    isValid = result["check"].as<bool>();
    list_size = result["list"].as<NODE_T>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  LightGraph G;
  Input input(inFname);
  input.readMtx(inFname,G);
  NODE_T n = G.numberOfNodes();
  //std::cout<<"graph reading complete"<<std::endl;
  //std::cout<<"n: "<< n <<" m: "<<G.numberOfEdges()<<std::endl;
  
  std::pair<NODE_T, NODE_T> maxD = getMaxDegreeNode(G);
  //std::cout<<"Maximum Degree: "<<maxD.second<<std::endl;
  if(list_size >=0)
    std::cout<<"Since list size is given, ignoring alpha"<<std::endl;
  PaletteColor palcol(n,target,alpha,list_size);

  
  double t1 = omp_get_wtime();
  for(NODE_T i =0 ;i <n;i++) {
    for(EDGE_T j=G.IA[i]; j < G.IA[i+1] ; j++) {
      if(i < G.JA[j]) {
        palcol.buildStreamConfGraph(i,G.JA[j]); 
      } 
    } 
  }  
  double createConfTime = omp_get_wtime() - t1;
  //std::cout<<"Conflict graph construction time: "<<createConfTime<<std::endl;


  std::vector< std::vector<NODE_T> > confEdges = palcol.getConfAdjList();
  //palcol.confColorGreedy();
  palcol.confColorLF();

  PalColStat palStat = palcol.getPalStat(); 
  palStat.m = G.numberOfEdges();
  palStat.confBuildTime = createConfTime;
  std::vector<NODE_T> colors = palcol.getColors();
  palStat.nColors = palcol.getNumColors(); 
  std::vector <NODE_T>  invVert = palcol.getInvVertices();
  printStat(0,palStat);
  std::cout<<"# of Invalid Vertices: "<<invVert.size()<<"\n";
  //std::cout<<"# of colors: " <<palcol.getNumColors()<<std::endl;
  if(isValid) {
    if(isValidColoring(G,colors))
      std::cout<<"Coloring valid"<<std::endl;
    else
      std::cout<<"Coloring invalid"<<std::endl;
  }
  return 0;
}  
