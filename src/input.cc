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


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>

#include "ClqPart/input.h"

#define DEBUG 1

bool static fileTypeCheck(std::string fn, std::string extension)
{
    if(fn.substr(fn.find_last_of(".") + 1) == extension) {
        return true;
    } else {
        return false;
    }
}

void Input::readGraphGen(std::string fileName, LightGraph &G, bool is_random,long seed, VAL_T minW, VAL_T maxW)
{
    if(fileTypeCheck(fileName,"mtx")==true)
    {
        readMtx(fileName,G,is_random,seed,minW,maxW);
    }
    else if(fileTypeCheck(fileName,"bin") == true)
    {
        readBinGen(fileName,G);
    }
    else
    {
        std::cout<<"file must be either mtx or bin."<<std::endl;
        std::exit(1);
    }
}


//read symmetric graph
void Input::readMtx(std::string fileName, LightGraph &G, bool is_random,long int seed, VAL_T minW, VAL_T maxW)
{
    //std::cout<<"start to read"<<std::endl;
    srand(seed);
    if(fileTypeCheck(fileName,"mtx")==false)
    {
        std::cout << "file type is not mtx"<<std::endl;
        std::exit(1);
    }
    std::ifstream fileread(fileName.c_str());
    if(fileread.is_open()==false)
    {
        std::cout << "No file named "<<fileName<<std::endl;
        std::exit(1);
    }

    //Read the mtx header to determine file type
    std::string header[6];
    fileread>>header[0]>>header[1]>>header[2]>>header[3]>>header[4];
    std::string weightType = header[3];
    fileread.seekg(0, std::ios::beg);


    //Ignore header and comments
    while (fileread.peek() == '%') fileread.ignore(2048, '\n');

    NODE_T nrow,ncol;
    EDGE_T nnz;
    fileread >> nrow >> ncol >> nnz;
    NODE_T n;
    n = nrow;

    std::vector<std::vector<NODE_T> > adjList;
    adjList.resize(n);

    NODE_T u;
    NODE_T v;
    VAL_T weight;
    EDGE_T nEdge = 0;

    for(EDGE_T i=0;i<nnz;i++) {
        if(weightType == "pattern") {
          fileread >> v >> u ;
        }
        else {
          fileread >> v >> u >> weight;
        }
        v--;
        u--;
        if(u<v) {
            adjList[u].push_back(v);
            //reverse edge
            adjList[v].push_back(u);
            nEdge++;
        }
    }
    G.setNumberEdges(nEdge);
    G.setNumberNodes(n);

    G.IA.resize(n+1);
    G.JA.resize(2*nEdge);

    //second pass. Build the three arrays
    G.IA[0] = 0;
    EDGE_T k=0;
    for(NODE_T i=0;i<n;i++)
    {
        for(NODE_T j=0;j<adjList[i].size();j++)
        {
            G.JA[k] = adjList[i][j];
            k++;
        }
        G.IA[i+1]= G.IA[i] + adjList[i].size();

        adjList[i].clear();
        std::vector<NODE_T>().swap(adjList[i]);

    }
    ///std::cout<<"exiting.."<<std::endl;

    //close the file
    fileread.close();

}


void Input::readBinGen(std::string fileName, LightGraph &G)
{
    std::ifstream inf;
    NODE_T nA,nB;
    EDGE_T m;

    if(fileTypeCheck(fileName,"bin")==false)
    {
        std::cout << "file type is not bin"<<std::endl;
        std::exit(1);
    }
    inf.open(fileName.c_str(),std::ios::in|std::ios::binary);
    if(inf.is_open()!=1) {
        std::cout<<"Can not open file "<<fileName<<std::endl;
        exit(1);                                
    }
    inf.read((char*)&nA, sizeof(NODE_T));
    inf.read((char*)&nB, sizeof(NODE_T));
    inf.read((char*)&m, sizeof(EDGE_T));
    //of.write((char*)&max, sizeof(unsigned int));

    G.IA.resize(nA+1);
    G.JA.resize(2*m);

    G.setNumberEdges(m);
    G.setNumberNodes(nA);

    inf.read((char*)&G.IA[0], sizeof(NODE_T) * (nA+1));
    inf.read((char*)&G.JA[0], sizeof(NODE_T) * 2*m);
    inf.close();

}

// EGR is the binary CSR format from Burtscher's ECL graph library (Texas
// State Univ.). Layout (all little-endian, host-byte-order ints):
//   int32 nodes
//   int32 edges            // length of nlist == 2 * (undirected edge count)
//   int32 nindex[nodes+1]  // CSR row pointers
//   int32 nlist[edges]     // CSR column indices (each undirected edge stored twice)
//   int32 eweight[edges]   // optional; we ignore weights for coloring
void Input::readEGR(std::string fileName, LightGraph &G)
{
    if(fileTypeCheck(fileName,"egr")==false)
    {
        std::cout << "file type is not egr"<<std::endl;
        std::exit(1);
    }
    FILE *f = std::fopen(fileName.c_str(), "rb");
    if(f == nullptr) {
        std::cout << "Can not open file "<<fileName<<std::endl;
        std::exit(1);
    }

    int32_t nodes = 0;
    int32_t edges = 0;
    if(std::fread(&nodes, sizeof(int32_t), 1, f) != 1 ||
       std::fread(&edges, sizeof(int32_t), 1, f) != 1) {
        std::cout << "Failed to read EGR header from "<<fileName<<std::endl;
        std::fclose(f);
        std::exit(1);
    }
    if(nodes < 1 || edges < 0) {
        std::cout << "Invalid EGR header: nodes="<<nodes<<" edges="<<edges<<std::endl;
        std::fclose(f);
        std::exit(1);
    }

    NODE_T n = static_cast<NODE_T>(nodes);
    EDGE_T mTotal = static_cast<EDGE_T>(edges);

    G.setNumberNodes(n);
    G.setNumberEdges(mTotal / 2);
    G.IA.resize(n + 1);
    G.JA.resize(mTotal);

    // nindex (n+1 int32s) -> IA (vector<EDGE_T>). Width may differ, so go via temp.
    std::vector<int32_t> tmpIA(n + 1);
    if(std::fread(tmpIA.data(), sizeof(int32_t), n + 1, f) !=
       static_cast<size_t>(n + 1)) {
        std::cout << "Failed to read EGR nindex from "<<fileName<<std::endl;
        std::fclose(f);
        std::exit(1);
    }
    for(NODE_T i = 0; i <= n; i++) G.IA[i] = static_cast<EDGE_T>(tmpIA[i]);

    // nlist (edges int32s) -> JA. JA is vector<NODE_T>; if NODE_T is int the
    // direct fread is safe, otherwise widen via temp buffer.
    if(sizeof(NODE_T) == sizeof(int32_t)) {
        if(std::fread(G.JA.data(), sizeof(int32_t), mTotal, f) !=
           static_cast<size_t>(mTotal)) {
            std::cout << "Failed to read EGR nlist from "<<fileName<<std::endl;
            std::fclose(f);
            std::exit(1);
        }
    } else {
        std::vector<int32_t> tmpJA(mTotal);
        if(std::fread(tmpJA.data(), sizeof(int32_t), mTotal, f) !=
           static_cast<size_t>(mTotal)) {
            std::cout << "Failed to read EGR nlist from "<<fileName<<std::endl;
            std::fclose(f);
            std::exit(1);
        }
        for(EDGE_T i = 0; i < mTotal; i++) G.JA[i] = static_cast<NODE_T>(tmpJA[i]);
    }
    // Optional eweight tail is ignored.
    std::fclose(f);
}

void Input::wrtBinGen(std::string fileName, LightGraph &G)
{
    std::ofstream of;
    NODE_T nA = G.numberOfNodes();
    EDGE_T m = G.numberOfEdges();

    if(fileTypeCheck(fileName,"bin")==false)
    {
        std::cout << "file type is not bin"<<std::endl;
        std::exit(1);
    }
    of.open(fileName.c_str(),std::ios::out|std::ios::binary);
    if(of.is_open()!=1) {
        std::cout<<"Can not open file "<<fileName<<std::endl;
        exit(1);                                
    }

    of.write((char*)&nA, sizeof(NODE_T));
    of.write((char*)&nA, sizeof(NODE_T));
    of.write((char*)&m, sizeof(EDGE_T));
    //of.write((char*)&max, sizeof(unsigned int));
    of.write((char*)&G.IA[0], sizeof(NODE_T) * (nA+1));
    of.write((char*)&G.JA[0], sizeof(NODE_T) *2*m);
    of.close();

}

