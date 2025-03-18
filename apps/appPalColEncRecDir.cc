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

#include <iostream>
#include <cstring>
#include<fstream>
#include <random>
#include <filesystem>


#include "ClqPart/JsonGraph.h"
#include "ClqPart/paletteCol.h"
#include "ClqPart/input.h"
#include "cxxopts/cxxopts.hpp"
#include "ClqPart/utility.h"
#include "ClqPart/greedyColor.h"
#include "ClqPart/MemUsage.h"

void printStat( int level, PalColStat &palStat) {
  
  std::cout<<"***********Level "<<level<<"*******"<<std::endl;
  std::cout<<"Num Nodes: "<<palStat.n<<"\n";
  std::cout<<"Num Edges: "<<palStat.m<<"\n";
  std::cout<<"Avg. Deg.: "<<(double) 2*palStat.m/palStat.n<<"\n";
  std::cout<<"Palette Size: "<<palStat.palSz<<"\n";
  std::cout<<"List Size: "<<palStat.lstSz<<"\n";
  std::cout<<"Num Conflict Edges: "<<palStat.mConf<<"\n";
  std::cout<<"Conflict to Edge (%): "<<(double)palStat.mConf/palStat.m*100<<"\n";
  std::cout<<"Num Colors: " <<palStat.nColors<<std::endl;
  std::cout<<"Assign Time: "<<palStat.assignTime<<"\n";
  std::cout<<"Conf. Build Time: "<<palStat.confBuildTime<<"\n";
  std::cout<<"Conf. Color Time: "<<palStat.confColorTime<<"\n";
  std::cout<<"\n";

}


class LOG {
public:
  int seed;
  int nLevels;
  NODE_T n;
  EDGE_T m; 
  std::string problem_name;
  std::string log_file;
  std::string order_conf;
  NODE_T tot_palsize;
  NODE_T tot_listsize;
  NODE_T ncolors;
  EDGE_T tot_conf_edges;
  EDGE_T max_conf_edges;
  double tot_build_time;
  double tot_conf_color_time;
  double tot_assign_time;
  double tot_invalid_time;
  double max_memory;
  std::fstream fwrite;

  LOG(int _seed, int _nLevels, NODE_T _n, EDGE_T _m, std::string _problem_name,
      std::string _log_file, std::string _order_conf, NODE_T _tot_palsize, NODE_T _tot_listsize, NODE_T _ncolors, 
      EDGE_T _tot_conf_edges, double _tot_build_time, double _tot_conf_time, 
      double _tot_assgn_time, double _tot_inv_time, double _mem) {
      seed = _seed; 
      nLevels = _nLevels;
      n = _n;
      m = _m;
      problem_name = _problem_name;
      log_file = _log_file;
      order_conf = _order_conf;
      tot_palsize = _tot_palsize;
      tot_listsize = _tot_listsize;
      ncolors = _ncolors;
      tot_conf_edges = _tot_conf_edges;
      tot_build_time = _tot_build_time;
      tot_conf_color_time = _tot_conf_time;
      tot_assign_time = _tot_assgn_time;
      tot_invalid_time = _tot_inv_time;
      max_memory = _mem;
      
      if(log_file.empty() == false) {
        bool fexist = false;
        if(std::filesystem::exists(log_file)) {
           fexist = true; 
        }
        fwrite.open(log_file,std::fstream::out | std::fstream::app);
        if(fwrite.is_open() == false) {
          std::cout<<"Could not open the result output file"<<std::endl; 
        }   
        if(fexist == false) 
          fwrite<<"problem,seed,nlevels,n,m,order_conf,tot_palsz,tot_lstsz,ncols,tot_conf_edges,max_conf_edges,tot_build_tm,tot_conf_tm,tot_assgn_tm,tot_inv_tm,mem"<<std::endl;
      }
  }
  LOG(std::string _log_file) {
      log_file = _log_file;
      if(log_file.empty() == false) {
        bool fexist = false;
        if(std::filesystem::exists(log_file)) {
           fexist = true; 
        }
        fwrite.open(log_file,std::fstream::out | std::fstream::app);
        if(fwrite.is_open() == false) {
          std::cout<<"Could not open the result output file"<<std::endl; 
        }   
        if(fexist == false) 
          fwrite<<"problem,seed,nlevels,n,m,order_conf,tot_palsz,tot_lstsz,ncols,tot_conf_edges,max_conf_edges,tot_build_tm,tot_conf_tm,tot_assgn_tm,tot_inv_tm,mem"<<std::endl;
      }
  }
  LOG()=default;
  void create_log_file(std::string _log_file) {
      log_file = _log_file;
      if(log_file.empty() == false) {
        bool fexist = false;
        if(std::filesystem::exists(log_file)) {
           fexist = true; 
        }
        fwrite.open(log_file,std::fstream::out | std::fstream::app);
        if(fwrite.is_open() == false) {
          std::cout<<"Could not open the result output file"<<std::endl; 
        }   
        if(fexist == false) 
          fwrite<<"problem,seed,nlevels,n,m,order_conf,tot_palsz,tot_lstsz,ncols,tot_conf_edges,max_conf_edges,tot_build_tm,tot_conf_tm,tot_assgn_tm,tot_inv_tm,mem"<<std::endl;
      }
  }
  void append() {
    fwrite<<problem_name<<","<<seed<<","<<nLevels<<","<<n<<","<<m<<","
              <<order_conf<<","<<tot_palsize<<","<<tot_listsize<<","<<ncolors<<","<<tot_conf_edges<<","<<max_conf_edges
              <<","<<tot_build_time<<","<<tot_conf_color_time<<","<<tot_assign_time
              <<","<<tot_invalid_time<<","<<max_memory<<std::endl; 
         
  }

};


int main(int argC, char *argV[]) {
  
  cxxopts::Options options("palcolEr", "read json pauli string files and color the graph using palette coloring algorithm (CPU only version)"); 
  options.add_options()
    ("in,infile", "json file containing the pauli strings", cxxopts::value<std::string>())
    ("prob,problem", "problem name", cxxopts::value<std::string>()->default_value(""))
    ("out,outfile", "json file containing the groups after coloring", cxxopts::value<std::string>()->default_value(""))
    ("res,result", "result log file ", cxxopts::value<std::string>()->default_value(""))
    ("t,target", "palette size", cxxopts::value<double>())
    ("a,alpha", "coefficient to log(n) for list size", cxxopts::value<float>()->default_value("1.0"))
    //("s,stream", "use streaming construction", cxxopts::value<bool>()->default_value("false"))
    ("l,list", "use explicit list size", cxxopts::value<NODE_T>()->default_value("-1"))
    ("inv,ninv", "number of invalid vertices tolerance", cxxopts::value<NODE_T>()->default_value("100"))
    ("o,order", "RANDOM, LIST",cxxopts::value<std::string>()->default_value("LIST"))
    ("c,check", "check validity of coloring", cxxopts::value<bool>()->default_value("false"))
    ("r,recurse", "use recursive coloring", cxxopts::value<bool>()->default_value("false"))
    ("sd,seed", "use seed", cxxopts::value<int>()->default_value("123"))
    ("h,help", "print usage")
    ;

  std::string inFname,outFname,orderName,resFileName,probName;
  int seed;
  double target1;
  NODE_T target,list_size,nInv;
  float alpha;
  bool isValid,isRec;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    probName = result["prob"].as<std::string>();
    outFname = result["outfile"].as<std::string>();
    resFileName = result["result"].as<std::string>();
    orderName = result["order"].as<std::string>();
    target1 = result["target"].as<double>();
    alpha = result["alpha"].as<float>();
    seed = result["seed"].as<int>();
    //isStream = result["stream"].as<bool>();
    isValid = result["check"].as<bool>();
    isRec = result["recurse"].as<bool>();
    list_size = result["list"].as<NODE_T>();
    nInv = result["ninv"].as<NODE_T>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }
  std::cout<<"The grouping output: "<<outFname<<std::endl;
  LOG log;
  
  if(resFileName.empty() == false) {
    log.create_log_file(resFileName); 
  }
  log.seed = seed;
  log.problem_name = probName;
  ClqPart::JsonGraph jsongraph(inFname, false, true); 
  NODE_T n = jsongraph.numOfData();
  log.n = n;
  double nextFrac;
  auto baseline = getPeakRSS();
  if(target1 < 1) {
    std::cout<<"Using target as node percentage"<<std::endl; 
    target = NODE_T(n*target1);
    nextFrac = target1;
  }
  else {
    target = NODE_T(target1); 
    nextFrac = 1.0/8.0;
  }
  if(list_size >=0)
    std::cout<<"Since list size is given, ignoring alpha"<<std::endl;
  PaletteColor palcol(n,target,alpha,list_size,seed);

  //std::cout<<getPeakRSS() - baseline<<std::endl;
  
  int level = 0;
  //level 0 
  //**************Conflict graph construction*****
  //double t1 = omp_get_wtime();
  /*if(isStream){
    ClqPart::Edge e;
    while(jsongraph.nextEdge(e)) {
       palcol.buildStreamConfGraph(e.u,e.v); 
    
    }
  }*/
  palcol.buildConfGraph<std::vector<uint32_t>>(jsongraph);
  //std::cout<<getPeakRSS() - baseline<<std::endl;
  
  //double createConfTime = omp_get_wtime() - t1;
  //std::cout<<"Conflict graph construction time: "<<createConfTime<<std::endl;

  //std::vector< std::vector<NODE_T> > confEdges = palcol.getConfAdjList();
  //palcol.confColorGreedy();
  if(orderName == "RANDOM"){
    std::cout<<"using random ordering for conflict coloring"<<"\n";
    palcol.confColorRand();
  }
  else {
    std::cout<<"using list greedy ordering for conflict coloring"<<"\n";
    palcol.confColorGreedy();
  }
  //std::vector<NODE_T> colors = palcol.getColors();
  //std::cout<<getPeakRSS() - baseline<<std::endl;
  
  std::vector <NODE_T>  invVert = palcol.getInvVertices();
  PalColStat palStat = palcol.getPalStat(level); 
  palStat.m = jsongraph.getNumEdge();
  //palStat.mConf = palcol.getNumConflicts();
  palStat.nColors = palcol.getNumColors(); 

  log.order_conf = orderName;
  log.m = palStat.m;
  log.tot_listsize += palStat.lstSz;
  log.tot_palsize += palStat.palSz;
  log.tot_assign_time += palStat.assignTime;
  log.tot_build_time += palStat.confBuildTime;
  log.tot_conf_color_time += palStat.confColorTime;
  log.tot_conf_edges += palStat.mConf;
  log.max_conf_edges = palStat.mConf;

  printStat(level,palStat);
  
  if (isRec == true) {
    while(invVert.size() > nInv) { 
      jsongraph.resetNumEdge();
      level++;
      if(invVert.empty() == false) {
        /*
        if(invVert.size() > 40000) alpha = 3;
        else if (invVert.size() > 20000) alpha = 2; 
        else if (invVert.size() > 5000) alpha = 1.5; 
        else alpha = 1;
        */
        //std::cout <<nextFrac<<std::endl;
        palcol.reInit(invVert,invVert.size()*nextFrac,alpha);
        //palcol.reInit(invVert,target,alpha);
        palcol.buildConfGraph<std::vector<uint32_t>>(jsongraph,invVert);
        if(orderName == "RANDOM"){
          palcol.confColorRand(invVert);
        }
        else {
          palcol.confColorGreedy(invVert);
        }
      }
      invVert = palcol.getInvVertices();
      palStat = palcol.getPalStat(level); 
      palStat.m = jsongraph.getNumEdge();
      //palStat.mConf = palcol.getNumConflicts();
      palStat.nColors = palcol.getNumColors(); 

      log.tot_listsize += palStat.lstSz;
      log.tot_palsize += palStat.palSz;
      log.tot_assign_time += palStat.assignTime;
      log.tot_build_time += palStat.confBuildTime;
      log.tot_conf_color_time += palStat.confColorTime;
      log.tot_conf_edges += palStat.mConf;
      if(palStat.mConf > log.max_conf_edges)
        log.max_conf_edges = palStat.mConf;

      printStat(level,palStat);
    }
  }
  log.nLevels = level+1;
  if(invVert.empty() == false) {
    std::cout<<"Final Num invalid Vert: "<<invVert.size()<<"\n";
    palcol.naiveGreedyColor<std::vector<uint32_t>>(invVert, jsongraph, palcol.getNumColors());
    palStat = palcol.getPalStat(level); 
    std::cout<<"Naive Color TIme: "<<palStat.invColorTime<<"\n";
    log.tot_invalid_time = palStat.invColorTime;
  }
  
  std::cout<<"# of Final colors: " <<palcol.getNumColors()<<std::endl;
  log.ncolors = palcol.getNumColors();
  auto mem = getPeakRSS()-baseline;
  log.max_memory = mem;
  std::cout<<"Total memory usage(KB): " << mem <<std::endl;
  
  if(isValid) {
    if(palcol.checkValidity(jsongraph)) 
      std::cout<<"Coloring valid"<<std::endl;
    else
      std::cout<<"Coloring invalid"<<std::endl;
  }

  if(outFname.empty() == false) {
    //Be advised that using this fuction is going to increase memory consumption, since we are recreating the data array. This is only needed if we are encoding the data.
    jsongraph.computeDataArray(true);
    std::vector<NODE_T> cols = palcol.getColors();
    json jsonGrp = jsongraph.createColGroup(cols, palcol.getNumColors()); 
    std::string jsonString = jsonGrp.dump(4);
    std::ofstream outFile(outFname);
    if(outFile.is_open()) {
      outFile << jsonString;
      outFile.close();
    } else {
      std::cerr<< outFname<<" can not be opened" <<std::endl; 
    }
  }
  if(log.fwrite.is_open() ) {
     log.append();  
  }
  
  return 0;
}  
