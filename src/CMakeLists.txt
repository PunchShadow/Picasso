set(HEADER_LIST
  "${PROJECT_SOURCE_DIR}/include/ClqPart/Types.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/JsonGraph.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/graphColor.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/graph.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/input.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/utility.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/greedyColor.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/paletteCol.h"
  "${PROJECT_SOURCE_DIR}/include/ClqPart/cuPaletteCol.cuh"
    )

# Make an automatic library - will be static or dynamic based on user setting
add_library(ClqPart
  JsonGraph.cc
  graph.cc
  input.cc
  greedyColor.cc
  paletteColor.cc
  utility.cc
  cuPaletteCol.cu
  ${HEADER_LIST}
  )

set_target_properties(ClqPart PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
  target_link_libraries(ClqPart PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(ClqPart PUBLIC
  "${PROJECT_SOURCE_DIR}/include"
  )
target_link_libraries(ClqPart PUBLIC ${CUDA_LIBRARIES})
target_include_directories(ClqPart PUBLIC ${CUDA_INCLUDE_DIRS})

#target_link_libraries(ClqPart PUBLIC ${ColPack_LIBS})
#target_link_directories(ClqPart PUBLIC ${ColPack_LIBDIR})
#target_include_directories(ClqPart PUBLIC ${ColPack_INCDIR})

# All users of this library will need at least C++11
target_compile_features(ClqPart PUBLIC cxx_std_11)
#target_compile_options(ra PRIVATE -O2 -fopenmp)
target_compile_options(ClqPart PRIVATE -O3)
