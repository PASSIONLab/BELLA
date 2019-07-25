cmake_minimum_required(VERSION 2.8)
project(gerbil)

# Set Default Build Type to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

# Compile Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive -w")

# Include Files
set(HEADER_FILES
        include/gerbil/Application.h
        include/gerbil/Bundle.h
        include/gerbil/CpuHasher.h
        include/gerbil/config.h
        include/gerbil/global.h
        include/gerbil/debug.h
        include/gerbil/FailureBuffer.h
        include/gerbil/FastFile.h
        include/gerbil/FastParser.h
        include/gerbil/FastReader.h
        include/gerbil/GpuHasher.h
        include/gerbil/KmcWriter.h
        include/gerbil/KMer.h
        include/gerbil/KmerCountingHashTable.h
        include/gerbil/KmerDistributer.h
        include/gerbil/KmerHasher.h
        include/gerbil/SequenceSplitter.h
        include/gerbil/SuperReader.h
        include/gerbil/SuperWriter.h
        include/gerbil/SyncQueue.h
        include/gerbil/TempFile.h
        include/gerbil/TempFileStatistic.h
        include/gerbil/ThreadBarrier.h
        include/gerbil/types.h
        include/cuda_ds/CountingHashTable.h
        include/cuda_ds/CountingHashTable/_CountingHashTable.h
        include/cuda_ds/CountingHashTable/_Key.h
        include/cuda_ds/CountingHashTable/_KeyValuePair.h
        )

# Source Files
set(SOURCE_FILES_CPP
        src/gerbil/Application.cpp
        src/gerbil/Bundle.cpp
        src/gerbil/FastFile.cpp
        src/gerbil/FastParser.cpp
        src/gerbil/FastReader.cpp
        src/gerbil/KmcWriter.cpp
        src/gerbil/KmerDistributor.cpp
        src/gerbil/SequenceSplitter.cpp
        src/gerbil/SuperReader.cpp
        src/gerbil/SuperWriter.cpp
        src/gerbil/TempFile.cpp
        src/gerbil/debug.cpp
        src/gerbil/gerbil.cpp
        src/gerbil/global.cpp
        )

set(SOURCE_FILES_CUDA
        src/cuda_ds/AddKernel.cu
        src/cuda_ds/CompressKernel.cu
        )

# import mandatory libraries
find_package(Boost REQUIRED COMPONENTS system thread filesystem regex)
find_package(BZip2 REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Threads)
find_package(CUDA)

# set include directory (-I)
include_directories(
        ${Boost_INCLUDE_DIRS}
        ${BZip2_INCLUDE_DIRS}
        ${ZLIB_INCLUDE_DIRS}
)

# gather up libraries
set(MY_LIBS
        ${Boost_LIBRARIES}
        ${BZIP2_LIBRARIES}
        ${ZLIB_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT})

# if CUDA is located at the system
if (CUDA_FOUND)
    add_definitions(-DGPU)                                                   # enable CUDA flag
    include_directories(${CUDA_INCLUDE_DIRS})
    set(MY_LIBS ${MY_LIBS} ${CUDA_LIBRARIES})
    cuda_compile(CUDA_OBJECTS ${SOURCE_FILES_CUDA})
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# build executable
add_executable(gerbil ${HEADER_FILES} ${SOURCE_FILES_CPP} ${CUDA_OBJECTS})

# build convert
add_executable(toFasta src/gerbil/toFasta.cpp)

# link against external libraries
target_link_libraries(gerbil ${MY_LIBS})

# install binary
install(TARGETS gerbil DESTINATION bin)

