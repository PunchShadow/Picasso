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
#include<vector>
#include<algorithm>

#include "ClqPart/Types.h"

struct Edge
{
  NODE_T u;
  NODE_T v;
};



//to hold list of edges.
typedef std::vector<Edge> EdgeList;


class LightGraph
{
  
  //number of nodes in the graph
  NODE_T n; 
  //number of edges in the graph
  EDGE_T m;
  
  public:  
  //default constructor
  LightGraph()
  {
    n=0;
    m=0;
  }
  //constructor
  LightGraph(NODE_T n1, EDGE_T m1)
  {
    n = n1;
    m = m1;

    //allocate memory for A,IA and JA
    //graph is recorded as symmetric matrix. thats why we need 2*m
    IA.resize(n+1);
    JA.resize(2*m);
  }
  //delete the graph
  ~LightGraph()
  {
    n = 0;
    m = 0;
    //if(RA != 0) delete []RA;
    IA.resize(0);
    JA.resize(0);
  }
  
  
  //Graph is always in CSR format
  std::vector<EDGE_T> IA;
  std::vector<NODE_T> JA;
  std::vector<EDGE_T> RA;
  std::vector<EDGE_T> invRA;

  void setNumberEdges(EDGE_T m1)
  {
    m = m1;
  }
  EDGE_T numberOfEdges()
  {
    return m;
  }
 

  void setNumberNodes(NODE_T n1)
  {
    n = n1;
  }
  NODE_T numberOfNodes()
  {
    return n;
  }
  
  //void createEdgeList(EdgeList &);
  //this create edgelist and if bool=1 it also create the G.RA array
  void createEdgeList(EdgeList &,bool=0);
  void printGraph(NODE_T = 0);
  void createRefArray(bool=false);
  void sortGraph();
  
};

/*
class BptGraph
{
  //number of nodes in the graph
  NODE_T nA;
  NODE_T nB;
  //number of edges in the graph
  EDGE_T m;
  
  public:  
  //default constructor
  BptGraph()
  {
    A= 0;
    IA = 0;
    JA = 0;
  }
  //constructor
  BptGraph(NODE_T n1, NODE_T n2, EDGE_T m1)
  {
    nA = n1;
    nB = n2;
    m = m1;

    //allocate memory for A,IA and JA
    A = new VAL_T[m+5];
    IA = new NODE_T[nA+5];
    JA = new NODE_T[m];
    //RA is a reference array. Given a undirected graph as an edgeList. it create the references needed to go back and forth to the three arrays and edgeList.
    //RA = new EDGE_T[m];
  }
  //delete the graph
  ~BptGraph()
  {
    nA = 0;
    nB = 0;
    m = 0;
    //if(A!=0) delete []A;
    //if(IA!=0) delete []IA;
    //if(JA!=0)delete []JA;
  }
  
  
  //Graph is always in CSR format
  VAL_T *A;
  NODE_T *IA;
  NODE_T *JA;
  //EDGE_T *RA;
  void setNumberEdges(EDGE_T m1)
  {
    m = m1;
  }
  EDGE_T numberOfEdges()
  {
    return m;
  }
 

  void setNumberNodes(NODE_T n1,NODE_T n2)
  {
    nA = n1;
    nB = n2;
  }
  NODE_T numberOfNodesA()
  {
    return nA;
  }

  NODE_T numberOfNodesB()
  {
    return nB;
  }
  void printGraph();
  void sortGraph();
};
*/
