#ifndef __COMMON__H_
#define __COMMON__H_

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>
#include <sys/stat.h>
#include <omp.h>

#define MAX_KMER_SIZE 32

#define MAX_FILE_PATH 1000
#define NO_GZIP

typedef long long int lld;
typedef long long unsigned llu;

#ifdef _OPENMP
#define MYTHREAD omp_get_thread_num()
#define THREADS omp_get_num_threads()
#define MAXTHREADS omp_get_max_threads()
#else
#define MYTHREAD 0
#define THREADS 1
#define MAXTHREADS 1
#endif


#define serial_printf printf
#define WARN printf

#define INT_CEIL(numerator, denominator) (((numerator) - 1) / (denominator) + 1);

static inline void * calloc_chk(size_t nmemb, size_t size)
{
    void *x = calloc(nmemb, size);
    return x;
}

static inline int fclose_track(FILE *f)
{
    int ret = fclose(f);
    return ret;
}


static FILE *fopen_chk(const char *path, const char *mode)
{
    FILE *f = fopen(path, mode);
    return f;
}


static off_t get_file_size(const char *fname)
{
    struct stat s;
    if (stat(fname, &s) != 0) {
       fprintf(stdout,"could not stat %s: %s\n", fname, strerror(errno));
        return 0;
    }
    return s.st_size;
}

static char *get_basename(char *bname, const char *path)
{
    const char *slash_pos = rindex(path, '/');
    assert(strlen(path) > 0);
    assert(strlen(path) < MAX_FILE_PATH);
    if (!slash_pos) {
        strcpy(bname, path);
    } else {
        assert(*slash_pos == '/');
        assert(strlen(slash_pos+1) > 0);
        strcpy(bname, slash_pos+1);
   	}
    return bname;
}


// returns 1 when it created the directory, 0 otherwise, -1 if there is an error
static int check_dir(const char *path) {
    if (0 != access(path, F_OK)) {
        if (ENOENT == errno) {
            // does not exist
            // note: we make the directory to be world writable, so others can delete it later if we
            // crash to avoid cluttering up memory
            if (0 != mkdir(path, 0777) && 0 != access(path, F_OK)) {
                fprintf(stderr, "Could not create the (missing) directory: %s (%s)", path, strerror(errno));
                return -1;
            }
        }
        if (ENOTDIR == errno) {
            // not a directory
            fprintf(stderr, "Expected %s was a directory!", path);
            return -1;
        }
    } else {
        return 0;
    }
    assert( access(path, F_OK) == 0 );
    return 1;
}


#ifndef MAX_RANKS_PER_DIR
#define MAX_RANKS_PER_DIR 1000
#endif

// replaces the given path with a rank based path, inserting a rank-based directory
// example:  get_rank_path("path/to/file_output_data.txt", rank) -> "path/to/per_thread/<rankdir>/<rank>/file_output_data.txt"
// of if rank == -1, "path/to/per_thread/file_output_data.txt"
static char * get_rank_path(char *buf, int rank) {
    int pathlen = strlen(buf);
    char newPath[MAX_FILE_PATH];
    char *lastslash = strrchr(buf, '/');
    int checkDirs = 0;
    int thisDir;
    char *lastdir = NULL;
    if (pathlen+25 >= MAX_FILE_PATH) {
        fprintf(stderr, "File path is too long (max: %d): %s\n", MAX_FILE_PATH, buf);
        return NULL;
    }
    if (lastslash) *lastslash = '\0';
    if (rank < 0) {
        if (lastslash) {
            snprintf(newPath, MAX_FILE_PATH, "%s/per_thread/%s", buf, lastslash+1);
            checkDirs = 1;
        } else {
            snprintf(newPath, MAX_FILE_PATH, "per_thread/%s", buf);
            checkDirs = 1;
        }
    } else {
        if (lastslash) {
            snprintf(newPath, MAX_FILE_PATH, "%s/per_thread/%08d/%08d/%s", buf, rank/MAX_RANKS_PER_DIR, rank, lastslash+1);
            checkDirs = 3;
        } else {
            snprintf(newPath, MAX_FILE_PATH, "per_thread/%08d/%08d/%s", rank/MAX_RANKS_PER_DIR, rank, buf);
            checkDirs = 3;
        }
    }
    strcpy(buf, newPath);
    while(checkDirs > 0) {
        strcpy(newPath, buf);
        thisDir = checkDirs;
        while(thisDir--) {
            lastdir = strrchr(newPath, '/');
            if (!lastdir) { fprintf(stderr, "What is happening here?!?!\n"); return NULL; }
            *lastdir = '\0';
        }
        check_dir(newPath);
        checkDirs--;
    }
    return buf;
}

#if defined (__cplusplus)
#include <string>
static std::string getRankPath(const char *path, int rank) {
    char buf[MAX_FILE_PATH];
    strcpy(buf, path);
    get_rank_path(buf, rank);
    return std::string(buf);
}
#endif

#endif
