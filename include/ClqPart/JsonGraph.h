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

//#include "ColPack_headers/ColPackHeaders.h"

#include "Types.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <omp.h>
#include "nlohmann/json.hpp"

//using json = nlohmann::ordered_json;
using json = nlohmann::json;


namespace ClqPart {

  struct Edge {
    NODE_T u;
    NODE_T v;
  };

  struct Graph {
    std::vector<Edge> edgeList;
    NODE_T n {0};
    NODE_T m {0};
  };
  
  class JsonGraph {
    double generateTime;
    double writeTime;
    json data;
    json dataAr;
    //json::iterator it,it1,beginIt;
    NODE_T u,v;
    EDGE_T numEdgeCom;
    NODE_T numDataPoints;
    public:
      JsonGraph() {}

      JsonGraph(std::string inFile,bool stream=false):inputFile(inFile) {    
        numEdgeCom = 0;
        std::ifstream f(inputFile);
        if (!f.is_open()) {
          std::cout<< "failed to open "<< inputFile<< "\n";
          exit(1);
        }
        data = json::parse(f);  
        f.close();

        //dataAr = nlohmann::json::array();
        dataAr = json::array();
        for (auto& el : data.items()) {
          json pair = json::array();
          pair.push_back(el.key());
          pair.push_back(el.value());
          dataAr.push_back(pair);
        }
        
        numDataPoints = dataAr.size();
        //std::cout<<numDataPoints<<std::endl;
        //beginIt = data.begin();
        if(stream == true) {
          //it = data.begin(); 
          //it1 = std::next(it,1);
          u = 0;
          v = 1;
        }
      }
      //void advance();
      void nextIndices();
      bool is_an_edge(std::string, std::string);
      bool is_an_edge(NODE_T, NODE_T);
      bool nextEdge(Edge &);
      void ReadJsonAdjacencyGraph(); 
      void ReadConstructWriteGraph(std::string fileName); 
      void writeGraphMtx(std::string fileName);
      double getGenTime() {return generateTime;}
      double getWriteTime() {return writeTime;}
      NODE_T numOfData() { return numDataPoints;}
      NODE_T getNumEdge() {return numEdgeCom;}
      void resetNumEdge() {numEdgeCom=0;}
      void printData();

    protected:
      std::string inputFile;
      Graph graph;
    
  };
}
