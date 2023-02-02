#include "ClqPart/graphColor.h"
#include "cxxopts/cxxopts.hpp"
#include <iostream>

int main(int argC, char *argV[]) {
  
  cxxopts::Options options("gcol", "read mtx and color the graph"); 
  options.add_options()
    ("in,infile", "input mtx file name", cxxopts::value<std::string>())
    ("out,outfile", "output color file name", cxxopts::value<std::string>())
    ("o,order", "LARGEST_FIRST,SMALLEST_LAST,NATURAL,RANDOM,DYNAMIC_LARGEST_FIRST,INCIDENCE_DEGREE",cxxopts::value<std::string>()->default_value("SMALLEST_LAST"))
    ("m,method","DISTANCE_ONE,DISTANCE_ONE_OMP",cxxopts::value<std::string>()->default_value("DISTANCE_ONE"))
    ("h,help", "print usage")
    ;

  std::string inFname, outFname,method,order;
  try{
    auto result = options.parse(argC,argV);
    if (result.count("help")) {
          std::cout<< options.help()<<"\n";
          std::exit(0);
    }
    inFname = result["infile"].as<std::string>();
    outFname = result["outfile"].as<std::string>();
    order = result["order"].as<std::string>(); 
    method = result["method"].as<std::string>(); 
  }
  catch(cxxopts::exceptions::exception &exp) {
    std::cout<<options.help()<<std::endl;
    exit(1);
  }
  
  colorGeneralGraph(inFname,outFname,order,method);

  return 0;
}
