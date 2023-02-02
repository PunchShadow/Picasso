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

#include "ClqPart/graphColor.h"
#include <cstring>
#include <unordered_set>

using namespace ColPack;

bool static fileTypeCheck(std::string fn, const std::string &extension) {
    if(fn.substr(fn.find_last_of(".") + 1) == extension) {
        return true;
    } else {
        return false;
    }
}

ColorStat colorGeneralGraph( std::string inFile, std::string outColorFile, std::string order, std::string method) {

  if(fileTypeCheck(inFile,"mtx")==false)
  {
    std::cout << inFile<<" input file type is not mtx"<<std::endl;
    std::exit(1);
  }

  GraphColoringInterface *g = new GraphColoringInterface(SRC_FILE, inFile.c_str(), "AUTO_DETECTED");
  g->Coloring(order.c_str(),method.c_str());

  double t1 = g->GetVertexOrderingTime();
  double t2 = g->GetVertexColoringTime();
  NODE_T nColor = g->GetVertexColorCount();
  ColorStat cstat{t1,t2,nColor};


  std::ofstream myfile(outColorFile.c_str());

  std::vector<int> vColors;
  g->GetVertexColors(vColors);

  if(myfile.is_open()) {
    for(NODE_T i=0;i<vColors.size();i++) { 
      myfile<< i <<" "<<vColors[i]<<"\n";
    }
    myfile.close();
  }
  else {
    std::cout<<"unable to open file named "<<outColorFile<<std::endl;
    exit(1);
  }

  return cstat;
}
