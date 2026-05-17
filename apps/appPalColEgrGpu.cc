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

// GPU palette colouring driver for general graphs in EGR (Burtscher CSR
// binary) format. Mirrors palcolEgr but funnels through the CUDA build
// path: PaletteColor<OffsetTy>::buildConfGraphGpuMemConscious(LightGraph)
// and the existing confColorGreedyCSR / confColorRandCSR.

#include <iostream>
#include <fstream>
#include <filesystem>
#include <climits>
#include <random>

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
  std::cout<<"Palette Size: "<<palStat.palSz<<"\n";
  std::cout<<"List Size: "<<palStat.lstSz<<"\n";
  std::cout<<"Num Conflict Edges: "<<palStat.mConf<<"\n";
  std::cout<<"Num Colors: "<<palStat.nColors<<std::endl;
  std::cout<<"Assign Time: "<<palStat.assignTime<<"\n";
  std::cout<<"Conf. Build Time: "<<palStat.confBuildTime<<"\n";
  std::cout<<"Conf. Color Time: "<<palStat.confColorTime<<"\n";
  std::cout<<"\n";
}

template<typename OffsetTy>
int runColoring(LightGraph &G, NODE_T target, float alpha, NODE_T list_size,
                int seed, double nextFrac, NODE_T nInv,
                const std::string &orderName, bool isRec, bool isValid,
                const std::string &outFname) {
  NODE_T n = G.numberOfNodes();
  PaletteColor<OffsetTy> palcol(n, target, alpha, list_size, seed);

  int level = 0;
  palcol.buildConfGraphGpuMemConscious(G);
  if(orderName == "RANDOM") palcol.confColorRandCSR();
  else                      palcol.confColorGreedyCSR();

  std::vector<NODE_T> invVert = palcol.getInvVertices();
  PalColStat palStat = palcol.getPalStat(level);
  palStat.nColors = palcol.getNumColors();
  printStat(level, palStat);

  if(isRec) {
    while(invVert.size() > (size_t)nInv) {
      level++;
      if(!invVert.empty()) {
        palcol.reInit(invVert, invVert.size() * nextFrac, alpha);
        palcol.buildConfGraphGpuMemConscious(G, invVert);
        if(orderName == "RANDOM") palcol.confColorRandCSR(invVert);
        else                      palcol.confColorGreedyCSR(invVert);
      }
      invVert = palcol.getInvVertices();
      palStat = palcol.getPalStat(level);
      palStat.nColors = palcol.getNumColors();
      printStat(level, palStat);
    }
  }

  if(!invVert.empty()) {
    std::cout<<"Final Num invalid Vert: "<<invVert.size()<<"\n";
    palcol.naiveGreedyColor(invVert, G, palcol.getNumColors());
    palStat = palcol.getPalStat(level);
    std::cout<<"Naive Color Time: "<<palStat.invColorTime<<"\n";
  }

  // Pure compute = assign(loop) + confBuild(kernel/CUB/scatter/sort) +
  // confColor + naive, summed over every level. copyTime (device alloc +
  // H2D/D2H memcpy) is reported separately and excluded.
  double pureCompute = 0.0, gpuCopy = 0.0;
  for(int lv = 0; lv <= level; lv++) {
    PalColStat s = palcol.getPalStat(lv);
    pureCompute += s.assignTime + s.confBuildTime + s.confColorTime + s.invColorTime;
    gpuCopy     += s.copyTime;
  }
  std::cout<<"Pure Compute Time: "<<pureCompute<<std::endl;
  std::cout<<"GPU Copy/Alloc Time: "<<gpuCopy<<std::endl;

  std::cout<<"# of Final colors: "<<palcol.getNumColors()<<std::endl;
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
  return palcol.getNumColors();
}

int main(int argC, char *argV[]) {
  cxxopts::Options options("palcolEgrG",
      "GPU palette coloring on EGR graphs (Burtscher CSR binary)");
  options.add_options()
    ("in,infile",  "EGR input file", cxxopts::value<std::string>())
    ("prob,problem","problem name",  cxxopts::value<std::string>()->default_value(""))
    ("out,outfile","plain-text output (vertex color per line)", cxxopts::value<std::string>()->default_value(""))
    ("res,result", "result log file (CSV, append)", cxxopts::value<std::string>()->default_value(""))
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
    if(result.count("help")) { std::cout<<options.help()<<"\n"; std::exit(0); }
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
  } catch(cxxopts::exceptions::exception &e) {
    std::cout<<options.help()<<std::endl;
    std::exit(1);
  }

  // Pay the one-time CUDA context init up front so it is NOT folded into
  // the first cudaMalloc inside assignListColor (which would pollute the
  // pure-compute measurement).
  double tw = omp_get_wtime();
  cudaFree(0);
  std::cout<<"CUDA Warmup Time: "<<(omp_get_wtime() - tw)<<std::endl;

  Input input;
  LightGraph G;
  double tl = omp_get_wtime();
  input.readEGR(inFname, G);
  std::cout<<"EGR Load Time: "<<(omp_get_wtime() - tl)<<std::endl;
  NODE_T n = G.numberOfNodes();
  EDGE_T m = G.numberOfEdges();
  std::cout<<"Loaded "<<inFname<<": n="<<n<<" m="<<m<<std::endl;

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
  if(target < 1) target = 1;
  if(list_size >= 0) std::cout<<"Since list size is given, ignoring alpha"<<std::endl;

  // Pick offset width: 32-bit covers up to ~2B conflict edges (each
  // undirected edge contributes 2 to atomicAdd targets and the COO write).
  // EGR's wire format is itself int32, so 2*m always fits in 32 bits.
  bool Offset32 = (uint64_t)m * 2 < (uint64_t)UINT_MAX;
  int finalColors = 0;
  if(Offset32) {
    std::cout<<"Using 32-bit offsets"<<std::endl;
    finalColors = runColoring<unsigned int>(G, target, alpha, list_size, seed,
                                            nextFrac, nInv, orderName, isRec,
                                            isValid, outFname);
  } else {
    std::cout<<"Using 64-bit offsets"<<std::endl;
    finalColors = runColoring<unsigned long long>(G, target, alpha, list_size, seed,
                                                  nextFrac, nInv, orderName, isRec,
                                                  isValid, outFname);
  }

  auto mem = getPeakRSS() - baseline;
  std::cout<<"Total memory usage(KB): "<<mem<<std::endl;

  if(!resFileName.empty()) {
    bool fexist = std::filesystem::exists(resFileName);
    std::ofstream resf(resFileName, std::ios::app);
    if(resf.is_open()) {
      if(!fexist) resf<<"problem,n,m,seed,palette,alpha,recurse,colors,mem_kb\n";
      std::string pName = probName.empty() ? getLastPartOfFilepath(inFname) : probName;
      resf<<pName<<","<<n<<","<<m<<","<<seed<<","<<target<<","<<alpha
          <<","<<(isRec?1:0)<<","<<finalColors<<","<<mem<<"\n";
      resf.close();
    }
  }
  return 0;
}
