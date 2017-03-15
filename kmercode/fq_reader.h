/* 
 * This can be used in either mpi or upc. Set USE_UPC_FOR_COMMON or USE_MPI_FOR_COMMON.
 */

#ifndef __FQ_READER_H
#define __FQ_READER_H

#include <stdint.h>

#ifndef NO_GZIP
#include "zlib.h"
#endif

#include "Buffer.h"

#ifndef FQ_READER_BUFFER_SIZE
#define FQ_READER_BUFFER_SIZE 8192
#endif

#if defined (__cplusplus)
extern "C" { 
#endif

struct fq_reader {
    FILE *f;
#ifndef NO_GZIP
    gzFile gz;
#endif
    int64_t fpos;
    int64_t size;
    int64_t start_read;
    int64_t end_read;
    // using long long instead of int64_t prevents cross-platform format errors with printf
    long long line;
    Buffer buf;
    char name[MAX_FILE_PATH];
    int max_read_len;
};

typedef struct fq_reader *fq_reader_t;

fq_reader_t create_fq_reader(void);
void destroy_fq_reader(fq_reader_t fqr);
void open_fq(fq_reader_t fqr, const char *fname, int cached_io);
void close_fq(fq_reader_t fqr);
//int get_next_fq_record_ptr(fq_reader_t fqr, char **id, char **nts, char **quals);
void hexifyId(char *name, int64_t *id1, int64_t *id2, int64_t step);
int get_next_fq_record(fq_reader_t fqr, Buffer id, Buffer nts, Buffer quals);
int get_fq_name_dirn(char *header, char **name, int *end);

int load_fq(fq_reader_t fqr, char *fname);
int unload_fq(char *fname);

int64_t estimate_fq(char *fname, int sampledReads, int64_t *estimatedReads, int64_t *estimatedBases);

#if defined (__cplusplus)
}
#endif

#endif
