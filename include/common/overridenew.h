#ifndef _OVERRIDE_NEW_
#define _OVERRIDE_NEW_

#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <tbb/scalable_allocator.h>

void* operator new(std::size_t sz) {
    return scalable_malloc(sz);
}
void operator delete(void* ptr) 
{
    scalable_free(ptr);
}

void* operator new[]( std::size_t sz ) {
    return scalable_malloc(sz);
}
void operator delete[]( void* ptr ) {
	scalable_free(ptr);
};

#endif
