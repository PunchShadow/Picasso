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

namespace ClqPart {

  struct Edge {
    NODE_T u;
    NODE_T v;
    VAL_T weight;
  };

  struct Graph {
    std::vector<Edge> edgeList;
    NODE_T n {0};
    NODE_T m {0};
  };
  
  class JsonGraph {
    public:
      JsonGraph() {}

      JsonGraph(std::string inFile)
        :inputFile(inFile) 
      {}

      void ReadJsonAdjacencyGraph(); 
      void writeGraphMtx(std::string fileName);

    protected:
      std::string inputFile;
      Graph graph;
    
  };
}
