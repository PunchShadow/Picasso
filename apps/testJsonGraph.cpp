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

//#include "nlohmann/json.hpp"
#include "ClqPart/JsonGraph.h"
#include "cxxopts/cxxopts.hpp"
#include <iostream>
#include <fstream>

//using json = nlohmann::json;

int main(int argC, char *argV[]) {
  cxxopts::Options options("testjson", "read json file"); 
  options.add_options()
    ("f,file", "file name", cxxopts::value<std::string>())
    ("h,help", "print usage")
    ;

  std::string fname;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    fname = result["file"].as<std::string>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }

  ClqPart::JsonGraph jsongraph(fname,true); 
  NODE_T n = jsongraph.numOfData();
  std::cout<<n<<std::endl;
  //jsongraph.ReadJsonAdjacencyGraph();
   
  return 0;
}
