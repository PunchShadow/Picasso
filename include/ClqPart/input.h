/*
 * Copyright (C) 2022  Ferdous,S M <ferdous.csebuet@egmail.com>
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

#include <string>

#include "ClqPart/graph.h"

class Input
{
  std::string fileName;

  public:
  void readMtx(std::string ,LightGraph&, bool=false, long int=0,VAL_T=0,VAL_T=0);
  void readBinGen(std::string , LightGraph &);
  void wrtBinGen(std::string , LightGraph &);
  void readGraphGen(std::string ,LightGraph&, bool=false, long int=0,VAL_T=0,VAL_T=0);

  Input()
  {
  
  }
  
  Input(std::string givenFile)
  {
    fileName = givenFile;
  }
};
