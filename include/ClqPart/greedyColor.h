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

class GreedyColor {
  std::string orderName;
  std::vector<NODE_T> vertexOrder; 
  std::vector<NODE_T> colors;
  LightGraph &G;
  NODE_T nColors;

public:

  GreedyColor(LightGraph &G1, std::string myorder = "LARGEST_FIRST"):G(G1) {
    NODE_T n = G.numberOfNodes();
    vertexOrder.resize(n);
    colors.resize(n,-1);
    nColors = 0;
    orderName = myorder;
  }

  void orderVertices(); 
  void colorVertices();
  NODE_T getNumColors() {return nColors;}
  std::vector<NODE_T> getColors() { return colors;}

private:
  void populateCandColors( NODE_T u, std::vector<NODE_T> &candColors); 
NODE_T firstAvailColor( NODE_T u, std::vector<NODE_T> &candColors)
  ;
  void orderLargestFirst();
};

