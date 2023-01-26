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


#include "ClqPart/JsonGraph.h"

#include "nlohmann/json.hpp"
#include <iostream>
#include <fstream>

using json = nlohmann::json;

namespace ClqPart {

  void JsonGraph::ReadJsonAdjacencyGraph() {

    std::ifstream f("/home/sferdou/ClqPartCpp/data/Ham1_H2O.json");
    json data = json::parse(f);  
    
    std::cout<<data;
  }
}
