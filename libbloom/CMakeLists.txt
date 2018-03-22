
message("Building libbloom")
add_subdirectory(murmur2)
include_directories(murmur2)

if(HIPMER_BLOOM64)
add_library(LIBBLOOM-NOUPC OBJECT bloom64.c)
add_library(LIBBLOOM OBJECT bloom64.c)
else()
add_library(LIBBLOOM-NOUPC OBJECT bloom.c)
add_library(LIBBLOOM OBJECT bloom.c)
add_library(LIBBLOOM OBJECT bloom.c)
endif()

if (UPC_ALLOCATOR)
     target_compile_definitions(LIBBLOOM PUBLIC -DUPC_ALLOCATOR=${UPC_ALLOCATOR})
endif()


