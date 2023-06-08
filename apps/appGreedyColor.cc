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

#include "ClqPart/greedyColor.h"
#include "ClqPart/input.h"
#include "cxxopts/cxxopts.hpp"


int main(int argC, char *argV[]) {
  
  cxxopts::Options options("greedycol", "read mtx and color the graph using greedy algorithm"); 
  options.add_options()
    ("in,infile", "input mtx file name", cxxopts::value<std::string>())
    ("out,outfile", "output color file name", cxxopts::value<std::string>()->default_value(""))
    ("o,order", "LARGEST_FIRST,SMALLEST_LAST,NATURAL,RANDOM,DYNAMIC_LARGEST_FIRST,INCIDENCE_DEGREE",cxxopts::value<std::string>()->default_value("LARGEST_FIRST"))
    ("h,help", "print usage")
    ;

  std::string inFname, outFname,order;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    outFname = result["outfile"].as<std::string>();
    order = result["order"].as<std::string>(); 
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  LightGraph G;
  Input input(inFname);
  input.readMtx(inFname,G);
  std::cout<<"graph reading complete"<<std::endl;
  std::cout<<"n: "<<G.numberOfNodes() <<" m: "<<G.numberOfEdges()<<std::endl;

  GreedyColor greedyCol(G);

  greedyCol.orderVertices();
  greedyCol.colorVertices();
  
  std::vector<NODE_T> colors = greedyCol.getColors();
  std::cout<<"coloring completed. num of colors: "<<greedyCol.getNumColors()<< " order + color time: "<< greedyCol.getOrderTime() <<" + "<<greedyCol.getColorTime()<<std::endl;

  if(outFname != "") {
    std::ofstream myfile(outFname);

    if(myfile.is_open()) {
      myfile<<greedyCol.getNumColors()<<"\n";
      for(NODE_T i=0;i<colors.size();i++) { 
        myfile<< i <<" "<<colors[i]<<"\n";
      }
      myfile.close();
    }
    else {
      std::cout<<"unable to open file named "<<outFname<<std::endl;
      exit(1);
    }
  }
  else {  
    std::cout<<greedyCol.getNumColors()<<"\n";
    for(NODE_T i=0;i<colors.size();i++) { 
      std::cout<< i <<" "<<colors[i]<<"\n";
    }
  }
  return 0;
}  
