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
    #ifdef ENABLE_GPU
    std::vector<uint32_t> dataEnc;
    #else
    std::vector<std::vector<uint32_t>> dataEnc;
    #endif
    size_t pauliEncSize;
    //json::iterator it,it1,beginIt;
    NODE_T u,v;
    EDGE_T numEdgeCom;
    NODE_T numDataPoints;
    public:
      JsonGraph() {}

      JsonGraph(std::string inFile,bool stream=false, bool encode=false):inputFile(inFile) {    
        numEdgeCom = 0;
        std::ifstream f(inputFile);
        if (!f.is_open()) {
          std::cout<< "failed to open "<< inputFile<< "\n";
          exit(1);
        }
        data = json::parse(f);  
        f.close();

        if(encode){
          constexpr uint32_t num_encoded_bits = sizeof(uint32_t) * 4;
          constexpr size_t terms_per_int = num_encoded_bits / 3;

          #ifdef ENABLE_GPU
          size_t pauli_length = data.begin().key().size();
          size_t num_strings = data.size();
          pauliEncSize = (pauli_length + terms_per_int - 1) / terms_per_int;
          dataEnc.reserve(pauli_length*pauliEncSize);
          #endif
          std::cout << "Encoding" << std::endl;
          for (auto& el : data.items()){
            // json pair = json::array();
            #ifndef ENABLE_GPU
            std::vector<uint32_t> matrix_encoding;
            #endif
            uint32_t cnt = 0;
            uint32_t encoded_int = 0;
            for (auto& matrix : el.key()) {
              if(cnt == num_encoded_bits/3){
                #ifdef ENABLE_GPU
                dataEnc.push_back(encoded_int);
                #else
                matrix_encoding.push_back(encoded_int);
                #endif
                cnt = 0;
                encoded_int = 0;
              }
              uint32_t mask = 0;
              // I = 000 bits, X = 110 bits, Y = 101 bits, Z = 011 bits
              if (matrix == 'I') {
                mask = 0b000;
              }
              else if (matrix == 'X') {
                mask = 0b110;
              }
              else if (matrix == 'Y') {
                mask = 0b101;
              }
              else if (matrix == 'Z') {
                mask = 0b011;
              }
              else {
                std::cout << "Invalid matrix: " << matrix << std::endl;
                exit(1);
              }
              encoded_int = encoded_int << 3;
              encoded_int = encoded_int | mask;
              cnt++;
            }
            #ifdef ENABLE_GPU
            dataEnc.push_back(encoded_int);
            #else
            matrix_encoding.push_back(encoded_int);
            dataEnc.push_back(matrix_encoding);
            #endif
          }
          numDataPoints = data.size();
          std::cout << "Done encoding!" << std::endl;
        }
        else{
          // dataAr = nlohmann::json::array();
          dataAr = json::array();
          for (auto& el : data.items()) {
            json pair = json::array();
            pair.push_back(el.key());
            pair.push_back(el.value());
            dataAr.push_back(pair);
          }
          numDataPoints = dataAr.size();
        }
        //beginIt = data.begin();
        if(stream == true) {
          //it = data.begin(); 
          //it1 = std::next(it,1);
          u = 0;
          v = 1;
        }
      }

      template <typename PauliTy = std::string>
      bool is_an_edge(NODE_T u, NODE_T v) {
        /*json::iterator bIt = data.begin();
        auto itU = std::next(bIt,u);
        auto itV = std::next(bIt,v);*/
        
        // if PauliTy is std::string
        if (std::is_same<PauliTy, std::string>::value) {
          PauliTy P = dataAr[u][0];
          PauliTy Q = dataAr[v][0];
          return is_an_edge(P,Q);
        }
        else{
          #ifdef ENABLE_GPU
          return is_an_edge(u,v);
          #else
          const std::vector<uint32_t> &P = dataEnc[u];
          const std::vector<uint32_t> &Q = dataEnc[v];
          return is_an_edge(P,Q);
          #endif
        }
      }

      template <typename PauliTy = std::string>
      bool nextEdge( Edge &e ) {
        while(1) {
          if(u >= numDataPoints-1) {
            return false;      
          }
          if (is_an_edge<PauliTy>(u,v)) {
            //advance();
            nextIndices();
            continue;
          }
          else {
            e.u = u;
            e.v = v;
            //advance();
            nextIndices();
            //numEdge++;
            //std::cout<<numEdge<<std::endl;
            return true;
          }

        }
        
      }

      //void advance();
      void nextIndices();
      bool is_an_edge(std::string, std::string);
      bool is_an_edge(const std::vector<uint32_t>&, const std::vector<uint32_t>&);
      #ifdef ENABLE_GPU
      bool is_an_edge(NODE_T P, NODE_T Q);
      #endif
      void ReadJsonAdjacencyGraph(); 
      void ReadConstructWriteGraph(std::string fileName); 
      void writeGraphMtx(std::string fileName);
      double getGenTime() {return generateTime;}
      double getWriteTime() {return writeTime;}
      NODE_T numOfData() { return numDataPoints;}
      NODE_T getNumEdge() {return numEdgeCom;}
      decltype(dataEnc)& getEncodedData() {return dataEnc;}
      size_t getPauliEncSize() {return pauliEncSize;}
      void printData();

    protected:
      std::string inputFile;
      Graph graph;
    
  };
}
