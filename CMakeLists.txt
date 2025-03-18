cmake_minimum_required(VERSION 3.9)
project(mps VERSION 0.1.0 LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # 最大最適化

# main executable
add_executable(${PROJECT_NAME}
  src/bucket.cpp
  src/loader.cpp
  src/main.cpp
  src/mps.cpp
  src/particle.cpp
  src/particles_exporter.cpp
  src/refvalues.cpp
  src/saver.cpp
  src/simulation.cpp
  src/weight.cpp
  src/divergence_calculator.cpp
  src/torque_calculator.cpp
  src/strain_calculator.cpp
  src/gradient_calculator.cpp
  src/laplacian_calculator.cpp
  src/neighbor_searcher.cpp
  src/Hamiltonian.cpp
)

# ------------------
# ----- OpenMP -----
# ------------------
find_package(OpenMP 4.5) #original OpenMP5.0 is required
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found. It runs on single thread.")
endif()

# -----------------
# ----- Eigen -----
# -----------------
include(FetchContent)
FetchContent_Declare(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
)
set(BUILD_TESTING OFF)
set(EIGEN_BUILD_TESTING OFF)
set(EIGEN_MPL2_ONLY ON)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(EIGEN_BUILD_DOC OFF)
FetchContent_MakeAvailable(Eigen)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)


# --------------------
# ----- yaml-cpp -----
# --------------------
include(FetchContent)
FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0
)
FetchContent_GetProperties(yaml-cpp)

if(NOT yaml-cpp_POPULATED)
  message(STATUS "Fetching yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()

target_link_libraries(${PROJECT_NAME} PUBLIC yaml-cpp::yaml-cpp)


# --------------------
# ----- argparse -----
# --------------------
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v3.0
)
FetchContent_MakeAvailable(argparse)
target_link_libraries(${PROJECT_NAME} PUBLIC argparse)
