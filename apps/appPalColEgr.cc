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

// Palette colouring driver for general graphs supplied in EGR (Burtscher's
// ECL graph) binary CSR format. Mirrors the recursive flow of palcolEr but
// skips the JsonGraph / Pauli encoding path entirely: the input file already
// contains the conflict graph, so we just iterate G's adjacency.

#include <iostream>
#include <cstring>
#include <fstream>
#include <random>
#include <filesystem>

#include "ClqPart/graph.h"
#include "ClqPart/paletteCol.h"
#include "ClqPart/input.h"
#include "cxxopts/cxxopts.hpp"
#include "ClqPart/utility.h"
#include "ClqPart/MemUsage.h"

void printStat(int level, PalColStat &palStat) {
  std::cout<<"***********Level "<<level<<"*******"<<std::endl;
  std::cout<<"Num Nodes: "<<palStat.n<<"\n";
  std::cout<<"Num Edges: "<<palStat.m<<"\n";
  std::cout<<"Avg. Deg.: "<<(double) 2*palStat.m/palStat.n<<"\n";
  std::cout<<"Palette Size: "<<palStat.palSz<<"\n";
  std::cout<<"List Size: "<<palStat.lstSz<<"\n";
  std::cout<<"Num Conflict Edges: "<<palStat.mConf<<"\n";
  if(palStat.m > 0)
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

  LOG()=default;
  void create_log_file(std::string _log_file) {
    log_file = _log_file;
    if(!log_file.empty()) {
      bool fexist = std::filesystem::exists(log_file);
      fwrite.open(log_file, std::fstream::out | std::fstream::app);
      if(!fwrite.is_open()) {
        std::cout<<"Could not open the result output file"<<std::endl;
      }
      if(!fexist)
        fwrite<<"problem,seed,nlevels,n,m,order_conf,tot_palsz,tot_lstsz,ncols,tot_conf_edges,max_conf_edges,tot_build_tm,tot_conf_tm,tot_assgn_tm,tot_inv_tm,mem"<<std::endl;
    }
  }
  void append() {
    fwrite<<problem_name<<","<<seed<<","<<nLevels<<","<<n<<","<<m<<","
          <<order_conf<<","<<tot_palsize<<","<<tot_listsize<<","<<ncolors<<","
          <<tot_conf_edges<<","<<max_conf_edges<<","<<tot_build_time<<","
          <<tot_conf_color_time<<","<<tot_assign_time<<","<<tot_invalid_time
          <<","<<max_memory<<std::endl;
  }
};

int main(int argC, char *argV[]) {
  cxxopts::Options options("palcolEgr",
      "Palette coloring on EGR-format graphs (Burtscher CSR binary; CPU)");
  options.add_options()
    ("in,infile",  "EGR input file", cxxopts::value<std::string>())
    ("prob,problem", "problem name", cxxopts::value<std::string>()->default_value(""))
    ("out,outfile","plain-text output (one 'vertex color' per line)", cxxopts::value<std::string>()->default_value(""))
    ("res,result", "result log file (CSV, append mode)", cxxopts::value<std::string>()->default_value(""))
    ("t,target",   "palette size: absolute (>=1) or fraction (0-1) of n", cxxopts::value<double>())
    ("a,alpha",    "coefficient to log(n) for list size", cxxopts::value<float>()->default_value("1.0"))
    ("l,list",     "explicit list size (overrides alpha if >=0)", cxxopts::value<NODE_T>()->default_value("-1"))
    ("inv,ninv",   "tolerance: stop recursion when invalid <= ninv", cxxopts::value<NODE_T>()->default_value("100"))
    ("o,order",    "RANDOM | LIST", cxxopts::value<std::string>()->default_value("LIST"))
    ("c,check",    "verify final coloring is proper", cxxopts::value<bool>()->default_value("false"))
    ("r,recurse",  "recursive refinement of invalid vertices", cxxopts::value<bool>()->default_value("false"))
    ("sd,seed",    "random seed", cxxopts::value<int>()->default_value("123"))
    ("h,help",     "print usage");

  std::string inFname, outFname, orderName, resFileName, probName;
  int seed;
  double target1;
  NODE_T target, list_size, nInv;
  float alpha;
  bool isValid, isRec;
  try {
    auto result = options.parse(argC, argV);
    if(result.count("help")) {
      std::cout<<options.help()<<"\n";
      std::exit(0);
    }
    inFname     = result["infile"].as<std::string>();
    probName    = result["prob"].as<std::string>();
    outFname    = result["outfile"].as<std::string>();
    resFileName = result["result"].as<std::string>();
    orderName   = result["order"].as<std::string>();
    target1     = result["target"].as<double>();
    alpha       = result["alpha"].as<float>();
    seed        = result["seed"].as<int>();
    isValid     = result["check"].as<bool>();
    isRec       = result["recurse"].as<bool>();
    list_size   = result["list"].as<NODE_T>();
    nInv        = result["ninv"].as<NODE_T>();
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    std::exit(1);
  }

  LOG log;
  if(!resFileName.empty()) log.create_log_file(resFileName);
  log.seed = seed;
  log.problem_name = probName.empty() ? getLastPartOfFilepath(inFname) : probName;

  Input input;
  LightGraph G;
  input.readEGR(inFname, G);
  NODE_T n = G.numberOfNodes();
  log.n = n;
  log.m = G.numberOfEdges();
  std::cout<<"Loaded "<<inFname<<": n="<<n<<" m="<<G.numberOfEdges()<<std::endl;

  auto baseline = getPeakRSS();

  double nextFrac;
  if(target1 < 1) {
    std::cout<<"Using target as node percentage"<<std::endl;
    target = NODE_T(n * target1);
    nextFrac = target1;
  } else {
    target = NODE_T(target1);
    nextFrac = 1.0/8.0;
  }
  if(target < 1) target = 1; // avoid zero-sized palette

  if(list_size >= 0)
    std::cout<<"Since list size is given, ignoring alpha"<<std::endl;

  PaletteColor palcol(n, target, alpha, list_size, seed);

  int level = 0;
  palcol.buildConfGraph(G);
  if(orderName == "RANDOM") {
    std::cout<<"using random ordering for conflict coloring"<<"\n";
    palcol.confColorRand();
  } else {
    std::cout<<"using list greedy ordering for conflict coloring"<<"\n";
    palcol.confColorGreedy();
  }

  std::vector<NODE_T> invVert = palcol.getInvVertices();
  PalColStat palStat = palcol.getPalStat(level);
  palStat.nColors = palcol.getNumColors();
  log.order_conf = orderName;
  log.tot_listsize    += palStat.lstSz;
  log.tot_palsize     += palStat.palSz;
  log.tot_assign_time += palStat.assignTime;
  log.tot_build_time  += palStat.confBuildTime;
  log.tot_conf_color_time += palStat.confColorTime;
  log.tot_conf_edges  += palStat.mConf;
  log.max_conf_edges  = palStat.mConf;
  printStat(level, palStat);

  if(isRec) {
    while(invVert.size() > (size_t)nInv) {
      level++;
      if(!invVert.empty()) {
        palcol.reInit(invVert, invVert.size() * nextFrac, alpha);
        palcol.buildConfGraph(G, invVert);
        if(orderName == "RANDOM") palcol.confColorRand(invVert);
        else                      palcol.confColorGreedy(invVert);
      }
      invVert = palcol.getInvVertices();
      palStat = palcol.getPalStat(level);
      palStat.nColors = palcol.getNumColors();
      log.tot_listsize    += palStat.lstSz;
      log.tot_palsize     += palStat.palSz;
      log.tot_assign_time += palStat.assignTime;
      log.tot_build_time  += palStat.confBuildTime;
      log.tot_conf_color_time += palStat.confColorTime;
      log.tot_conf_edges  += palStat.mConf;
      if(palStat.mConf > log.max_conf_edges) log.max_conf_edges = palStat.mConf;
      printStat(level, palStat);
    }
  }
  log.nLevels = level + 1;

  if(!invVert.empty()) {
    std::cout<<"Final Num invalid Vert: "<<invVert.size()<<"\n";
    palcol.naiveGreedyColor(invVert, G, palcol.getNumColors());
    palStat = palcol.getPalStat(level);
    std::cout<<"Naive Color Time: "<<palStat.invColorTime<<"\n";
    log.tot_invalid_time = palStat.invColorTime;
  }

  std::cout<<"# of Final colors: "<<palcol.getNumColors()<<std::endl;
  log.ncolors = palcol.getNumColors();
  auto mem = getPeakRSS() - baseline;
  log.max_memory = mem;
  std::cout<<"Total memory usage(KB): "<<mem<<std::endl;

  if(isValid) {
    if(palcol.checkValidity(G)) std::cout<<"Coloring valid"<<std::endl;
    else                        std::cout<<"Coloring invalid"<<std::endl;
  }

  if(!outFname.empty()) {
    std::ofstream outFile(outFname);
    if(outFile.is_open()) {
      std::vector<NODE_T> cols = palcol.getColors();
      for(NODE_T i = 0; i < n; i++) outFile<<i<<" "<<cols[i]<<"\n";
      outFile.close();
    } else {
      std::cerr<<outFname<<" can not be opened"<<std::endl;
    }
  }

  if(log.fwrite.is_open()) log.append();
  return 0;
}
