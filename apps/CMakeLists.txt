# Make an automatic library - will be static or dynamic based on user setting
#add_executable(cliqePart cliqePart.cpp)
#target_link_libraries(cliqePart PUBLIC ${ColPack_LIBS})
#target_link_directories(cliqePart PUBLIC ${ColPack_LIBDIR})
#target_include_directories(cliqePart PUBLIC ${ColPack_INCDIR})
#target_include_directories(cliqePart PUBLIC "${PROJECT_SOURCE_DIR}/include") 

#if(OpenMP_CXX_FOUND)
#  target_link_libraries(cliqePart PUBLIC OpenMP::OpenMP_CXX)
#endif()


add_executable(testjson testJsonGraph.cpp)
target_include_directories(testjson PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(testjson PUBLIC ClqPart)

add_executable(jsonmtx jsonToMtx.cpp)
target_include_directories(jsonmtx PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(jsonmtx PUBLIC ClqPart)

#add_executable(gcol appGraphColor.cc)
#target_include_directories(gcol PUBLIC "${PROJECT_SOURCE_DIR}/include") 
#target_link_libraries(gcol PUBLIC ClqPart)

add_executable(grdcol appGreedyColor.cc)
target_include_directories(grdcol PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(grdcol PUBLIC ClqPart)

add_executable(palcol appPaletteColor.cc)
target_include_directories(palcol PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(palcol PUBLIC ClqPart)

add_executable(palcolS appPalColDir.cc)
target_include_directories(palcolS PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(palcolS PUBLIC ClqPart)

add_executable(naiveCol appNaiveColDir.cc)
target_include_directories(naiveCol PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(naiveCol PUBLIC ClqPart)

add_executable(tGraph testGraphIn.cc)
target_include_directories(tGraph PUBLIC "${PROJECT_SOURCE_DIR}/include") 
target_link_libraries(tGraph PUBLIC ClqPart)

#message(STATUS "ColPack: found=${ColPack_FOUND}, include_dirs=${ColPack_INCDIR}, "
#  "lib=${ColPack_LIBDIR}, version=${ColPack_VERSION}")
  #target_include_directories(gstream PUBLIC
  #"${PROJECT_SOURCE_DIR}/include"
  #)

# All users of this library will need at least C++11
#target_compile_features(gstream PUBLIC cxx_std_17)
#target_compile_options(ra PRIVATE -O2 -fopenmp)
#target_compile_options(gstream PRIVATE -O3)
