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

#ifndef TYPES_H_
#define TYPES_H_

#include <iostream>
#include <cstdio>
#include "debug.h"
#include "config.h"

namespace gerbil {

typedef int int32;
typedef long int64;
typedef short int16;
typedef char int8;

typedef unsigned char uchar;
typedef uchar byte;				// Alias for uchar
typedef uchar uint8;			// Alias for uchar
typedef unsigned short ushort;
typedef ushort uint16;			// Alias for ushort
typedef unsigned int uint;
typedef uint uint32;			// Alias for uint
typedef unsigned long ulong;
typedef ulong uint64;			// Alias for ulong

typedef uint_fast16_t uint_tfn;	// uint for tempFilesNumber
typedef uint32 uint_cv;			// uint for counterValue
typedef uint_fast8_t uint_tid;	// uint for threadIds

#define TEMPFILEID_NONE static_cast<uint_fast16_t>(0xffff)

typedef enum {
	st_genome, st_reads
} TSeqType;

typedef enum {
	ft_unknown, ft_fasta, ft_fastq, ft_multiline
} TFileType;

typedef enum {
	fc_none, fc_gzip, fc_bz2, fc_DECOMPRESSOR
} TFileCompr;

typedef enum {
	ht_free, ht_locked, ht_occupied
} THashtableLockState;

typedef enum {
	of_none, of_gerbil, of_fasta
} TOutputFormat;

#define C_0  453569
#define C_1  5696063
#define C_2  3947602847
#define C_3  342971
#define C_4  6127577

//#define MEMORY_BARRIER asm volatile("" ::: "memory");
inline void memoryBarrier() {
	asm volatile("" ::: "memory");
}

#define CAS(addr, old, val) __sync_val_compare_and_swap(addr, old, val)

#if __BYTE_ORDER == __LITTLE_ENDIAN
//supported
#elif __BYTE_ORDER == __BIG_ENDIAN
#error unsupported byte order
#else
#error unknown byte order
#endif

}

#endif /* TYPES_H_ */
