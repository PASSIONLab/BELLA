/*********************************************************************************
Copyright (c) 2016 Marius Erbert, Steffen Rechner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*********************************************************************************/

#ifndef CONFIG_H_
#define CONFIG_H_

#include <assert.h>
#include <climits>

namespace gerbil {


#define VERSION_MAJOR 1
#define VERSION_MINOR 11

//#define DEB
//#define DEV
#define REL

#define HYBRID_COUNTER false

#ifdef GPU
#define IF_GPU(x) x
#else
#define IF_GPU(x)
#endif


#ifdef REL
#undef DEB
#undef DEV
#define IF_REL(x) x
#endif

#ifdef DEV
#undef DEB
#define IF_DEV(x) x
#define IF_DEB_DEV(x) x
#endif

#ifdef DEB
#define IF_DEB(x) x
#define IF_DEB_DEV(x) x
#endif

#ifndef REL
#define IF_REL(x)
#endif

#ifndef DEV
#define IF_DEV(x)
#endif

#ifndef DEB
#define IF_DEB(x)
#endif

#if !(defined(DEB) || defined(DEV))
#define IF_DEB_DEV(x)
#endif


#define FILL 0.4
#define FILL_MAX 0.8
#define FILL_GPU 0.5
#define START_RATIO FILL

#define HISTOGRAM_SIZE 512


/*
 * default parameters
 */
#define DEF_MINIMIZER_SIZE 7
#define MIN_MINIMIZER_SIZE 5
#define MAX_MINIMIZER_SIZE 9

#define DEF_KMER_RANGE 128
#define DEF_KMER_SIZE 28
#define MIN_KMER_SIZE 8
#define MAX_KMER_SIZE DEF_KMER_RANGE+MIN_KMER_SIZE

#define DEF_MEMORY_SIZE ((uint64_t)  4 * 1024)
#define MIN_MEMORY_SIZE ((uint64_t)       512)
#define MAX_MEMORY_SIZE ((uint64_t) 1024 * 1024 * 1024)

#define DEF_TEMPFILES_NUMBER 512
#define MIN_TEMPFILES_NUMBER 2
#define MAX_TEMPFILES_NUMBER 4 * 1024

#define DEF_THREADS_NUMBER 8
#define MIN_THREADS_NUMBER 4
#define MAX_THREADS_NUMBER 128

#define DEF_THRESHOLD_MIN 3
#define MIN_THRESHOLD_MIN 1
#define MAX_THRESHOLD_MIN (1 << 30)

#define DEF_NORM true
#define DEF_VERBOSE IF_DEB_DEV(true)IF_REL(false)


#define FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD 1

#define SB_WRITER_THREADS_NUMBER 1


#define GB_TO_B(x) ((uint64_t)x << 30)
#define MB_TO_B(x) ((uint64_t)x << 20)
#define KB_TO_B(x) ((uint64_t)x << 10)

#define B_TO_GB(x) ((uint64_t)x >> 30)
#define B_TO_MB(x) ((uint64_t)x >> 20)
#define B_TO_KB(x) ((uint64_t)x >> 10)

// sizes of bundles
#define FAST_BLOCK_SIZE_B                KB_TO_B( 64)
#define FAST_BUNDLE_DATA_SIZE_B        KB_TO_B(512)
#define READ_BUNDLE_SIZE_B                KB_TO_B(128)
#define SUPER_BUNDLE_DATA_SIZE_B        KB_TO_B( 32)
#define KMER_BUNDLE_DATA_SIZE_B            KB_TO_B(128)
#define KMC_BUNDLE_DATA_SIZE_B            KB_TO_B(256)
#define GPU_KMER_BUNDLE_DATA_SIZE        KB_TO_B(256)

// step1: size of memory
#define RUN1_MEMORY_GENERAL_B            MB_TO_B( 32)
#define MIN_FASTBUNDLEBUFFER_SIZE_B        MB_TO_B( 64)
#define MIN_READBUNDLEBUFFER_SIZE_B        MB_TO_B( 64)
#define MIN_SUPERBUNDLEBUFFER_SIZE_B    MB_TO_B( 64)
#define MIN_SUPERWRITERBUFFER_SIZE_B    MB_TO_B( 64)
#define MAX_FASTBUNDLEBUFFER_SIZE_B        MB_TO_B(128)
#define MAX_READBUNDLEBUFFER_SIZE_B        MB_TO_B(256)
#define MAX_SUPERBUNDLEBUFFER_SIZE_B    MB_TO_B(256)

// step2: size of memory
#define RUN2_MEMORY_GENERAL_B            MB_TO_B( 48)
#define MIN_SUPERBUNDLEBUFFER2_SIZE_B    MB_TO_B( 32)
#define MIN_KMERBUNDLEBUFFER_SIZE_B    MB_TO_B( 32)
#define MIN_KMCHASHTABLE_SIZE_B        MB_TO_B( 32)
#define MIN_KMCBUNDLEBUFFER_SIZE_B        MB_TO_B( 32)
#define GPU_COPY_BUFFER_SIZE            MB_TO_B( 16)

#define MEM_KEY_HT    0.8

#define NULL_BUCKET_VALUE UINT_MAX


#define LOOP2(N, X) X(N - 1); X(N)
#define LOOP4(N, X) LOOP2(N-2, X); LOOP2(N, X)
#define LOOP8(N, X) LOOP4(N-4, X); LOOP4(N, X)
#define LOOP16(N, X) LOOP8(N-8, X); LOOP8(N, X)
#define LOOP32(N, X) LOOP16(N-16, X); LOOP16(N, X)
#define LOOP64(N, X) LOOP32(N-32, X); LOOP32(N, X)
#define LOOP128(N, X) LOOP64(N-64, X); LOOP64(N, X)
#define LOOP256(N, X) LOOP128(N-128, X); LOOP128(N, X)
#define LOOP512(N, X) LOOP256(N-256, X); LOOP256(N, X)

}

#endif /* CONFIG_H_ */
