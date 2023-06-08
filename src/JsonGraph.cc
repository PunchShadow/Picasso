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
#include <cmath>
#include <complex>
#include <string>

using json = nlohmann::ordered_json;

namespace ClqPart {
  std::pair<std::complex<double>,std::string> pqMerge(std::string P, std::string Q)
  {
    auto L = P.size();
    EDGE_T nimag = 0; // numberf of 'i'
    EDGE_T nsign = 0; // number of '-1'
    std::string new_string = "";
    for (auto i = 0; i < L; i++)
    {
        if (P[i] == 'I')
        {
            new_string += Q[i];
        }
        else if (Q[i] == 'I')
        {
            new_string += P[i];
        }
        else if (P[i] == Q[i])
        {
            new_string += 'I';
        }
        else if (P[i] == 'X')
        {
            if (Q[i] == 'Y')
            {
                new_string += 'Z';
                nimag += 1;
            }
            else if (Q[i] == 'Z')
            {
                new_string += 'Y';
                nimag += 1;
                nsign += 1;
            }
        }
        else if (P[i] == 'Y')
        {
            if (Q[i] == 'X')
            {
                new_string += 'Z';
                nimag += 1;
                nsign += 1;
            }
            else if (Q[i] == 'Z')
            {
                new_string += 'X';
                nimag += 1;
            }
        }
        else if (P[i] == 'Z')
        {
            if (Q[i] == 'X')
            {
                new_string += 'Y';
                nimag += 1;
            }
            else if (Q[i] == 'Y')
            {
                new_string += 'X';
                nimag += 1;
                nsign += 1;
            }
        }
    }
    //
    std::complex<double> sgn = std::pow(std::complex<double>(0, 1), nimag) * std::pow(-1, nsign);
    //
    return std::make_pair(sgn,new_string);
  }

  bool is_an_edge(std::string P, std::string Q) {
    int L = P.length();
    int cnt = 0;
    for (int i = 0; i < L; i++) {
      if (P[i] != 'I' && Q[i] != 'I' && P[i] != Q[i]) {
        cnt++;
      }
    }
    if (cnt % 2 == 1) {
      return true;
    }
    else {
      return false;
    }
  }

  bool static fileTypeCheck(std::string fn, std::string extension)
{
    if(fn.substr(fn.find_last_of(".") + 1) == extension) {
        return true;
    } else {
        return false;
    }
}

  void JsonGraph::ReadJsonAdjacencyGraph() {

    std::ifstream f(inputFile);
    if (!f.is_open()) {
      std::cout<< "failed to open "<< inputFile<< "\n";
      exit(1);
    }
    json data = json::parse(f);  
    //json::iterator it = data.begin(); 
    
    /*
    for (auto& [key, val] : data.items())
    {
      std::string s(val);
      std::cout << std::fixed<< "key: " << key << ", value:" << std::stod(s) << '\n';
    } 
    */
    
    double t1 = omp_get_wtime();     
    auto i=0;
    std::cout<<graph.edgeList.max_size()<<std::endl;
    for (auto it = data.begin(); it != std::prev(data.end(),1); ++it, ++i)
    {
      auto j=i+1;
      for (auto it1 = std::next(it,1); it1 != data.end(); ++it1,++j) {
        //std::string s(it.value());
        //auto v1 = std::complex<double>(std::stod(s),0);
        //std::string s1(it1.value());
        //auto v2 = std::complex<double>(std::stod(s1),0);
        //std::cout << "key: " << it.key() << ", value:" << std::stod(s) << '\n';
        //std::cout << "key: " << it1.key() << ", value:" << std::stod(s1) << '\n';

        //std::cout<<"\n";
        //auto Z1 = pqMerge(it.key(),it1.key());
        //auto Z2 = pqMerge(it1.key(),it.key());

        //auto _tmp = Z1.first*std::conj(v1)*v2+Z2.first*std::conj(v2)*v1;

        /*if(Z1.second == Z2.second && std::abs(_tmp.real())<1e-6) {
          continue;
        }*/
        if (is_an_edge(it.key(),it1.key())) {
            continue;
        }
        else {
          //std::cout<<Z1.second<<" "<<Z2.second<<" "<<i+1<<" "<<j+1<<"\n";
          Edge e{i,j};
          graph.edgeList.push_back(e);
          if(graph.edgeList.size() % 10000000 == 0) {
            std::cout<<graph.edgeList.size()<<std::endl; 
          }
        }
      }
    }
    graph.n = i+1;
    graph.m = graph.edgeList.size();
    generateTime = omp_get_wtime() - t1;
  }


  void JsonGraph::writeGraphMtx(std::string fileName) {
    if(fileTypeCheck(fileName,"mtx")==false)
    {
        std::cout << "file type is not mtx"<<std::endl;
        std::exit(1);
    }
    std::ofstream myfile(fileName.c_str());

    double t1 = omp_get_wtime();     
    if(myfile.is_open())
    {
        myfile << "%%MatrixMarket matrix coordinate pattern symmetric"<<std::endl;
        myfile << graph.n<<" "<<graph.n<<" "<<graph.m<<"\n";
        for(auto elem: graph.edgeList) { 
          myfile<< elem.v+1<<" "<<elem.u+1<<"\n";
        }
        myfile.close();
    }
    else
        std::cout<<"unable to open file"<<std::endl;

    writeTime = omp_get_wtime() - t1;
     
  }
  void JsonGraph::ReadConstructWriteGraph(std::string fileName) {

    std::ifstream f(inputFile);
    if (!f.is_open()) {
      std::cout<< "failed to open "<< inputFile<< "\n";
      exit(1);
    }
    if(fileTypeCheck(fileName,"mtx")==false) {
      std::cout << "file type is not mtx"<<std::endl;
      std::exit(1);
    }

    std::ofstream myfile(fileName.c_str());

    if(myfile.is_open() == false) {
      std::cout<<"unable to open file"<<std::endl; 
      exit(1);
    }
    json data = json::parse(f);  
    //json::iterator it = data.begin(); 
    
    /*
    for (auto& [key, val] : data.items())
    {
      std::string s(val);
      std::cout << std::fixed<< "key: " << key << ", value:" << std::stod(s) << '\n';
    } 
    */
    
    double t1 = omp_get_wtime();     
    NODE_T i=0;
    EDGE_T m = 0;
    for (auto it = data.begin(); it != std::prev(data.end(),1); ++it, ++i)
    {
      auto j=i+1;
      for (auto it1 = std::next(it,1); it1 != data.end(); ++it1,++j) {
        if (is_an_edge(it.key(),it1.key())) {
            continue;
        }
        else {
          //std::cout<<Z1.second<<" "<<Z2.second<<" "<<i+1<<" "<<j+1<<"\n";
          m++;
          myfile<<j+1<<" "<<i+1<<"\n";
        }
      }
    }
    NODE_T n = i+1;
    myfile<<n<<" "<<n<<" "<<m<<std::endl;
    generateTime = omp_get_wtime() - t1;
  }
}
