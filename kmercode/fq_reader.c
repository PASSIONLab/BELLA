/*
 * Note: this fastq reader code is not specific to upc. It is also used by kMerCount with mpi.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <stdint.h>
#include <errno.h>
#include <sys/time.h>

#include "fq_reader.h"

#ifndef MMAP_BLOCK_SIZE
#define MMAP_BLOCK_SIZE getpagesize()
#endif

#ifndef DEFAULT_READ_LEN
#define DEFAULT_READ_LEN 100000
#endif

#ifdef __MACH__
#include <mach/mach_time.h>

static double InitConversionFactor() {
  mach_timebase_info_data_t timebase;
  mach_timebase_info(&timebase);
  double conversion_factor = (double)timebase.numer / (double)timebase.denom;
  return conversion_factor;
}

static double now() {
  static double conversion_factor = 0.0;
  if (conversion_factor == 0.0) conversion_factor = InitConversionFactor();
  return mach_absolute_time() * conversion_factor;
}
 
#else

static double now() {
   struct timespec ts;
   clock_gettime(CLOCK_REALTIME, &ts);
   return ts.tv_sec + ts.tv_nsec / 1000000.0;
}

#endif

fq_reader_t create_fq_reader(void)
{ 
    fq_reader_t fqr = (fq_reader_t) calloc_chk(1, sizeof(struct fq_reader));
    fqr->buf = initBuffer(DEFAULT_READ_LEN);
    assert(isValidBuffer(fqr->buf));
    return fqr;
}

void destroy_fq_reader(fq_reader_t fqr)
{
    if (fqr) {
        if (fqr->buf) 
            freeBuffer(fqr->buf);
        fqr->buf = NULL;

#ifndef NO_GZIP
        if (fqr->gz)
            close_fq(fqr);
        fqr->gz = NULL;
#endif
        if (fqr->f) 
            close_fq(fqr);
        fqr->f = NULL;
        fqr->size = 0;
        free(fqr);
    }
}


static void strip_trailing_spaces(char *s)
{
    char *end = s + strlen(s) - 1;
    while (end >= s && isspace(*end))
        end--;
    *(end + 1) = '\0';
}

inline static char *get_fq_name(char *header)
{
    assert(header != NULL);
    assert(header[0] == '@');
    int len = strlen(header);
//    if (len >= MAX_READ_NAME_LEN - 1) {
//       return NULL;
//    }
    strip_trailing_spaces(header);
    len = strlen(header);
    // only convert if new illumina 2 format
    if (header[len - 2] != '/') {
        // Latest Illumina header format
        char *end = strchr(header, '\t');
        if (!end) {
            end = strchr(header, ' ');
            if (!end) {
                // no comment, just return the name without modification
                return header;
            }
        }
        // truncate name at comment
        assert(*end == ' ' || *end == '\t');
        end[0] = '\0';
        // check for @pair/1 @/pair2 vs @pair 1:Y:... @pair 2:Y:....
        int commentpos = end - header;
        if (commentpos > 3 && header[commentpos-2] == '/' && (header[commentpos-1] == '1' || header[commentpos-1] == '2')) {
            // @pair/1 or @pair/2
            return header;
        } else if ((len < commentpos + 7) || end[2] != ':' || end[4] != ':' || end[6] != ':' || (end[1] != '1' && end[1] != '2')) {
            // unknown pairing format
            return header;
        }
        // @pair 1:Y:.:.: or @pair 2:Y:.:.:
        // replace with @pair/1 or @pair/2
        end[0] = '/';
        end[2] = '\0';
    }
    assert(header[0] == '@');
    return header;
}

int get_fq_name_dirn(char *header, char **name, int *end)
{
    *name = get_fq_name(header);
    if (!*name) {
        return 0;
    }
    int len = strlen(*name);
    if (len < 3) {
        return 0;
    }
    char end_char = (*name)[len - 1];
    if ((*name)[len-2] != '/' || (end_char != '1' && end_char != '2')) {
        return 0;
    }
    *end = (end_char == '1' ? 1 : 2);
    (*name)[len - 2] = '\0';
    return 1;
}




// returns 0 for unpaired (offset to '@'), or the offset to '1' or '2' for paired
static int8_t get_pairidx_from_name_line(char *s)
{
    // Supports two types of paired formats:
    //
    // Illumina:
    // @name/1
    // @name/2  
    // 
    // Casava 1.8:
    // @HISEQ15:136:C6G2YANXX:4:1101:1074:2037 1:N:0:AGTTCC
    // @HISEQ15:136:C6G2YANXX:4:1101:1074:2037 2:N:0:AGTTCC
    // 
    // Further addition:
    // We've seen data in the Casava format (the second one above) without the final sequence 
    // fragment. So we support that also now.
    
    assert(s[0] == '@');
    int pairIdx = 0, max = 2047;
    int endIdx = 1, spaceIdx = 0;
    while (s[endIdx] != '\n' && s[endIdx] != '\0') {
        if (endIdx >= max) 
            break;
        if (spaceIdx == 0 && isspace(s[endIdx])) 
            spaceIdx = endIdx;
        endIdx++;
    }
    // line too long to determine pairing!
    if (s[endIdx] != '\n') 
        return 0; 
    // something is wrong - the line is too short
    if (endIdx <= 3)
        return 0; 
    // Illumina with random comment
    if (spaceIdx > 3 && s[spaceIdx - 2] == '/' && (s[spaceIdx - 1] == '1' || s[spaceIdx - 1] == '2'))
        return spaceIdx - 1;
    // Illumina without comment
    if (s[endIdx - 2] == '/' && (s[endIdx - 1] == '1' || s[endIdx - 1] == '2'))
        return endIdx - 1;
    // Casava
    if (s[spaceIdx + 2] == ':' && s[spaceIdx + 4] == ':' && s[spaceIdx + 6] == ':' && 
        (s[spaceIdx + 1] == '1' || s[spaceIdx + 1] == '2') &&
        (s[spaceIdx + 3] == 'Y' || s[spaceIdx + 3] == 'N')) 
        return spaceIdx + 1;
    return 0;
}

static int64_t get_fptr_for_next_record(fq_reader_t fqr, int64_t offset)
{
    if (offset == 0) return 0; // first record is the first record, include it.  Every other partition will be at least 1 full record after offset.
#ifndef NO_GZIP
    if(fqr->gz || !fqr->f) printf(stderr,"Can not fseek in a compressed file! (%s)\n", fqr->name);
#else
    if(!fqr->f) printf(stderr,"Must fseek in an open file! (%s)\n", fqr->name);
#endif
    char lastName[256];
    lastName[0] = '\0';
    if (offset >= fqr->size) return fqr->size;
    int ret = fseek(fqr->f, offset, SEEK_SET);
    if (ret != 0) {
        printf(stderr,"fseek could not execute on %s to %lld: %s", fqr->name, (lld) offset, strerror(errno));
        return -1;
    }

    printf(stdout,"Finding fptr for next record of %s after %lld\n", fqr->name, (lld) offset);

    // skip first (likely parial) line after this offset to ensure we start at the beginning of a line
    if (!fgetsBuffer(fqr->buf, 2047, fqr->f)) 
        return ftell(fqr->f);
    int64_t new_offset = 0;
    int count = 0, last_pair_line = 0;
    char prev_pair = 0;
    int64_t prev_offset = ftell(fqr->f);
    Buffer readLines = initBuffer(4096);
    printfBuffer(readLines, "%s", getStartBuffer(fqr->buf));
    int possibleHeaderLine = -1;
    char prevName[256]; prevName[0] = '\0';
    while (1) {
        new_offset = ftell(fqr->f);
        resetBuffer(fqr->buf);
        if (!fgetsBuffer(fqr->buf, 2047, fqr->f))
            break;
        count++;
        char *buf = getStartBuffer(fqr->buf);
        printfBuffer(readLines, "%s", buf);
        //fprintf(stdout, "Read line of %lld length at offset %lld count %d: %s", (lld) getLengthBuffer(fqr->buf), (lld) prev_offset, count, buf);
        if (buf[0] == '@') {
            // '@' immediately after newline yields two possiblities: This is a header line or a quality line (or the file is invalid).  Use logic to determine 
            // if it is a quality line, the next line must be a header
            // if it is a header line the 4 lines from now will be a header too (3 lines may also be a quality line with starting '@' too)
            // once the first header is determined, look for interleaved pairs, but allow unpaired reads
            if (possibleHeaderLine > 0 && count == possibleHeaderLine + 1) {
               // possibleHeaderLine was actually a quality line.  so this MUST be a header line
               possibleHeaderLine = count; // the next header line we expect will be 4 after (and should never come back here)
               //LOGF("header line after quality line at offset=%lld count=%d: %s", (lld) new_offset, count, buf);
               prev_pair = 0; // reset any bogus pair
               prevName[0] = '\0'; // reset the previous name
            } else if (possibleHeaderLine > 0 && count == possibleHeaderLine + 4) {
               // the first possibleHeaderLine was actualy a header.  good from here. this MUST be a header line
               possibleHeaderLine = count; // the next header line we expect will be 4 after (and will hit here again)
               //LOGF("header line after 4 lines after first header at offset=%lld count=%d: %s", (lld) new_offset, count, buf);
            } else if (possibleHeaderLine > 0 && count == possibleHeaderLine + 3) {
               // this must be a quality line after a real header line
               //LOGF("Found quality line after 3 lines after the first header offset=%lld count=%d: %s", (lld) new_offset, count, buf);
               continue; // not interested in this
            } else {
               // first '@' line this could be either header or quality lets skip it
               if (possibleHeaderLine > 0) fprintf(stderr,"This should not happen!  possibleHeaderLine=%d at %lld.  Read from %lld offset: '%s'\n", possibleHeaderLine, (lld) new_offset, (lld) offset, getStartBuffer(readLines));
               possibleHeaderLine = count;
               //LOGF("possibleHeaderLine may be quality line at offset=%lld count=%d: %s\n", (lld) new_offset, count, buf);
               // Try to use this line anyway, even if it happens to be quality
            }
            // now look to pairs (or not)
            strncpy(lastName, buf, 255);
            int pairIdx = get_pairidx_from_name_line(buf);
            char pair = 'X'; // 'X' for unpaired, '\0' for not assigned yet, '1', '2' for read 1, read 2 respectively
            if (pairIdx) {
              pair = buf[pairIdx];
            }
            //LOGF("Discovered pair: '%c' lastName=%s prevName=%s\n", pair, lastName, prevName);
            if (!prev_pair) {
                prev_pair = pair;
                prev_offset = new_offset;
                last_pair_line = count;
                //DBG("no prev_par at count %d\n", count);
            } else {
                if (last_pair_line + 1 == count) {
                    // last pair was actually a quality line
                    prev_pair = pair;
                    prev_offset = new_offset;
                    last_pair_line = count;
                    //LOGF("Last line was actually a quality line using next at %lld: '%s'\n", (lld) new_offset, buf);
                    fprintf(stderr, "This should no longer happen! %s\n", getStartBuffer(readLines));
                } else if (last_pair_line + 4 == count) {
                    // not interleaved or interleaved and prev pair was 1 and this is 2
                    if (prev_pair == '1' && pair == '2' && strcmp(prevName, lastName) == 0)  {
                        new_offset = prev_offset;
                        //LOGF("Interleaved this record completes the pair, using prev_offset %lld at count %d, lastName: %s prevName: %s\n", (lld) prev_offset, count, lastName, prevName);
                    } else if (prev_pair == '2' && pair == '1' && strcmp(prevName, lastName) != 0) {
                        //LOGF("Interleaved this record starts a new pair, using new_offset %lld at count %d, lastName: %s prevName: %s\n", (lld) new_offset, count, lastName, prevName);
                    } else {
                        new_offset = prev_offset;
                        //LOGF("Not interleaved at count %d, using prev_offset %lld, lastName: %s prevName: %s\n", count, (lld) prev_offset, lastName, prevName);
                    }
                    break;
                }
            }
            strncpy(prevName, lastName, 255);
        }
        if (count > 13) 
            fprintf(stderr,"Could not find a valid line in the fastq file %s, last line: %s\n", fqr->name, buf);
    }
    assert(new_offset >= offset);
    //fprintf(stdout,"Found record at %lld (count=%d %lld bytes after requested offset=%lld):\n%s...Skipped: %.*s\n", (lld) new_offset, count, (lld) (new_offset - offset), (lld) offset, getStartBuffer(readLines) + (new_offset - offset), (int) (new_offset - offset), getStartBuffer(readLines));
    freeBuffer(readLines);
    // for interleawed, make sure to the next record starts with a /1
    //fprintf(stdout,"Chose %lld as my offset '%s' (prev_offset: %lld, startOffset: %lld)\n", (lld) new_offset, lastName, (lld) prev_offset, (lld) offset);
    return new_offset;
}

void open_fq(fq_reader_t fqr, const char *fname, int cached_io) {
    fqr->line = 0;
    fqr->max_read_len = 0;
    // we have a single file for all threads
    strcpy(fqr->name, fname);
    fqr->size = get_file_size(fqr->name);
 

    fqr->f = fopen_chk(fqr->name, "r");
    int64_t read_block = INT_CEIL(fqr->size, THREADS);
    fqr->start_read = read_block * MYTHREAD;
    fqr->end_read = read_block * (MYTHREAD + 1);
    if (MYTHREAD > 0) 
	    fqr->start_read = get_fptr_for_next_record(fqr, fqr->start_read);
    if (MYTHREAD == THREADS - 1)
            fqr->end_read = fqr->size;
    else 
            fqr->end_read = get_fptr_for_next_record(fqr, fqr->end_read);

    //long long int start_read = fqr->start_read;
    //long long int end_read = fqr->end_read;    	
    //fprintf(stdout,"thread %d of %d and max of %d start block %lld and end block %lld out of file size %lld, each read block is %lld \n",MYTHREAD, THREADS, MAXTHREADS, start_read, end_read, fqr->size, read_block );
    
    
    assert(fqr->end_read >= fqr->start_read);
    if (setvbuf(fqr->f, NULL, _IOFBF, FQ_READER_BUFFER_SIZE) != 0) 
        fprintf(stdout,"Could not setvbuf on %s to %d! %s", fname, FQ_READER_BUFFER_SIZE, strerror(errno));
    if (fseek(fqr->f, fqr->start_read, SEEK_SET) != 0) {
        long long int start_read = fqr->start_read;
        fprintf(stderr,"reset_my_partition could not fseek on %s to %lld: %s\n", fqr->name, start_read, strerror(errno));
    }
    fqr->fpos = fqr->start_read;
    assert(fqr->fpos == ftell(fqr->f));
    if(MYTHREAD==0)
    	fprintf(stdout, "reading FASTQ file %s\n", fname);
#ifndef __APPLE__
    if (fqr->f) {
        posix_fadvise(fileno(fqr->f), fqr->start_read, fqr->end_read - fqr->start_read, 
                      POSIX_FADV_SEQUENTIAL);
    }
#endif
}

// read just a portion of the file (from fqr->start to fqr->end) and store to /dev/shm
// this function assums fqr is already at the start position
int load_fq(fq_reader_t fqr, char *fname)
{
    open_fq(fqr, fname, 0); // never cached_io so the file is read in partitions
    assert(fqr->f);
    int64_t my_fsize = fqr->end_read - fqr->start_read;
    //fprintf(stdout,"load_fq storing %lld bytes from %s (%lld - %lld)\n", (lld) my_fsize, fqr->name, (lld) fqr->start_read, (lld) fqr->end_read);
    int err = 0;
    // create shm file
    char my_fname[MAX_FILE_PATH];
    char bname[MAX_FILE_PATH];
    sprintf(my_fname, "/dev/shm/%s", get_rank_path(get_basename(bname, fname), MYTHREAD));
#ifndef NO_GZIP
    sprintf(my_fname, "%s.gz", my_fname);
    int blockSize = 8*1024*1024;
    char *buf = (char*) malloc_chk(blockSize);
    GZIP_FILE f = gzopen_chk(my_fname, (YES_GZIP ? "w1" : "w"));
    int64_t bytes_remaining = my_fsize;
    double t1, t2;
    double readTime = 0.0, writeTime = 0.0;
    t1 = now();
    do {
        int64_t readBlock = bytes_remaining < blockSize ? bytes_remaining : blockSize;
        int64_t bytesRead = fread(buf, 1, readBlock, fqr->f);
        t2 = now();
        readTime += t2-t1;
        //fprintf(stdout,"load_fq read %lld bytes from %s\n", (lld) bytesRead, fqr->name);
        int64_t bytesWritten = 0;
        do {
           int64_t thisWrite = gzwrite(f, buf+bytesWritten, bytesRead - bytesWritten);
           //fprintf(stdout,"load_fq wrote %lld uncompressed bytes to %s\n", (lld) thisWrite, my_fname);
           if (!thisWrite) fprintf(stderr,"Could not write %lld bytes to %s\n", (lld) bytesRead - bytesWritten, my_fname);
           bytesWritten += thisWrite;
        } while (bytesWritten < bytesRead);
        bytes_remaining -= bytesRead;
        t1 = now();
        writeTime += t1 - t2;
    } while (bytes_remaining > 0);
    gzclose(f);
    t2 = now();
    writeTime += t2-t1;
    free(buf);
    int64_t compressed_size = get_file_size(my_fname);
    
    // store the uncompressed size in a secondary file
    sprintf(my_fname, "%s.uncompressedSize", my_fname);
    FILE *fsizefd = fopen_chk(my_fname, "w");
    fwrite(&my_fsize, sizeof(int64_t), 1, fsizefd);
    fclose_track(fsizefd);
    //fprintf(stdout,"Compressed %lld bytes to %lld (to %0.3f %%) in %0.3f read %0.3f write\n", (lld) my_fsize, (lld) compressed_size, 100.0 * compressed_size / my_fsize, readTime, writeTime);

#else // def NO_GZIP
    int fd = open(my_fname, O_CREAT|O_TRUNC|O_RDWR, S_IRUSR|S_IWUSR);
    if (fd != -1) {
        if (ftruncate(fd, my_fsize) == -1) {
            fprintf(stdout,"Failed to truncate %s: %s\n", my_fname, strerror(errno));
            err = errno;
        } else if (my_fsize > 0) {
            //fprintf(stdout,"Opened shm file %s for writing of size %lld\n", my_fname, (lld) my_fsize);
            void *tmpaddr = mmap(NULL, my_fsize, PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
            if (tmpaddr == MAP_FAILED) {
                fprintf(stdout,"Could not mmap %s for writing: %s\n", my_fname, strerror(errno));
                err = errno;
            } else {
                size_t bytes_read = 0, bytes_remaining = my_fsize;
                while (bytes_remaining > 0) {
                    int64_t bytes = fread((char*)tmpaddr + bytes_read, 1, bytes_remaining, fqr->f);
                    if (bytes == 0) 
                        break;
                    bytes_read += bytes;
                    bytes_remaining -= bytes;
                }
                if (bytes_read != my_fsize) {
                    fprintf(stdout,"Read the wrong number of bytes: %lld vs %ld for %s\n",
                         (lld) bytes_read, my_fsize, my_fname); 
                    err = 1;
                }
                if (msync(tmpaddr, my_fsize, MS_SYNC) == -1) {
                    fprintf(stdout,"msync returned error %s\n", strerror(errno));
                    err = errno;
                }
                if (munmap(tmpaddr, my_fsize) == -1) {
                    fprintf(stdout,"munmap returned error %s\n", strerror(errno));
                    err = errno;
                }
            }
        }
        if (close(fd) != 0)
            fprintf(stdout,"Could not close the shm file %s! %s", my_fname, strerror(errno));
        //else
            //fprintf(stdout,"Wrote %s of %lld size.\n", my_fname, (lld) my_fsize);
    } else {
        err = errno;
    }
#endif // NO_GZIP
    close_fq(fqr);
    return -err;
}

int unload_fq(char *fname)
{
    char my_fname[MAX_FILE_PATH];
    char bname[MAX_FILE_PATH];
#ifdef NO_GZIP
    sprintf(my_fname, "/dev/shm/%s", get_rank_path(get_basename(bname, fname), MYTHREAD));
#else
    sprintf(my_fname, "/dev/shm/%s.gz", get_rank_path(get_basename(bname, fname), MYTHREAD));
#endif
    if (unlink(my_fname) == -1) 
        return -errno;
    return 0;
}

// include newline char
static char *get_next_line(fq_reader_t fqr)
{
    char *buf = NULL;
    if (fqr->fpos > fqr->end_read) return NULL;
    size_t len = 0, maxLen = (fqr->end_read - fqr->fpos);
    assert(isValidBuffer(fqr->buf));
    resetBuffer(fqr->buf);
    assert(getStartBuffer(fqr->buf) == getEndBuffer(fqr->buf));
    if (fqr->fpos >= fqr->end_read) return NULL;
#ifndef NO_GZIP
    if (fqr->gz) {
       assert(! fqr->f);
       buf = gzgetsBuffer(fqr->buf, (maxLen > 2048 ? 2048 : maxLen), fqr->gz);
    } else
#endif
    buf = fgetsBuffer(fqr->buf, (maxLen > 2048 ? 2048 : maxLen), fqr->f);
    if (!buf)
        return NULL;
        
    len = getLengthBuffer(fqr->buf);
    fqr->fpos += len;
#ifndef NO_GZIP
    if (fqr->gz) {
        assert(fqr->fpos == gztell(fqr->gz));
    } else
#endif
    assert(fqr->fpos == ftell(fqr->f));
    assert(buf != NULL);
    assert(buf == getStartBuffer(fqr->buf));
    if (len == 0) return NULL;
    assert(buf[len] == '\0');
    if (buf[len-1] != '\n') {
        WARN("should be end of line, not '%c' %d: '%s' at line %lld\n", 
             buf[len], buf[len], buf, fqr->line);
    }
    fqr->line++;
    assert(len > 0);
    assert(len == getLengthBuffer(fqr->buf));
    return buf;
}

void hexifyId(char *name, int64_t *id1, int64_t *id2, int64_t step) {
    assert(name[0] == '@');
    int len = strlen(name);
    int pairlen = len;
    int pair = 0;
    int64_t *id = id1;
    if (len > 3 && name[len-2] == '/') {
       // has a pair
       pairlen -= 2;
       if (name[len-1] == '1') {
          pair = 1;
       } else if (name[len-1] == '2') {
          pair = 2;
          id = id2;
       }
    }
    char newname[32];
    sprintf(newname, "%llx%s", (long long int) *id, (pairlen==len?"":name + pairlen));
    strcpy(name+1, newname);
    *id += step;
}

int get_next_fq_record(fq_reader_t fqr, Buffer id, Buffer nts, Buffer quals) {
    resetBuffer(id);
    resetBuffer(nts);
    resetBuffer(quals);
    if (fqr->fpos >= fqr->end_read) {
        return 0;
    }
    Buffer lastBuf = fqr->buf;

    // get all four lines, one for each field
    for (int i = 0; i < 4; i++) {
        char *tmp = NULL;
        switch(i) {
            case 0: fqr->buf = lastBuf; break;
            case 1: fqr->buf = nts;     break; 
            case 2: fqr->buf = lastBuf; break;
            case 3: fqr->buf = quals;   break;
            default: assert(0);
        };
        tmp = get_next_line(fqr);
        if (!tmp) {
            // return the proper buffer back
            fqr->buf = lastBuf;
            return 0;
        }

        char *buf = getStartBuffer(fqr->buf);
        if (getLengthBuffer(fqr->buf) == 0 || buf[0] == '\n') {
            // empty
            if (i != 0 ) 
                fprintf(stderr,"Invalid FASTQ at line %lld, expected empty or just newline at the end: '%s'\n", (lld) fqr->line, buf);
            // return the proper buffer back
            fqr->buf = lastBuf;
            return 0;
        }
            
        if (i == 0) {
            // header line
            if (buf[0] != '@')        
                fprintf(stderr,"Invalid FASTQ at line %lld, expected read name (@):\n%s\n", (lld) fqr->line, buf);
            // construct universally formatted name (illumina 1 format)
            strip_trailing_spaces(buf);
            char *read_name = get_fq_name(buf);
            if (!read_name)
                fprintf(stderr,"Invalid FASTQ name format: %s\n", buf);
            strcpyBuffer(id, read_name);
        } else if (i == 1) {
            // sequence
            chompBuffer(fqr->buf);
#ifdef CONFIG_CHECK_SEQS        
            char label[256];
            sprintf(label, "  in %s:%d at line %lld\n", __FILE__, __LINE__, (lld) fqr->line);
            assert(strlen(buf) == getLengthBuffer(fqr->buf));
            int len = check_seqs(buf, label);
#endif
            assert(fqr->buf == nts);
        } else if (i == 2) {
            if (buf[0] != '+')
                fprintf(stderr,"Invalid FASTQ at line %lld, expected '+':\n%s\n", (lld) fqr->line, buf);
        } else if (i == 3) {
            chompBuffer(fqr->buf);
            assert(fqr->buf == quals);
        } else {
            fprintf(stderr,"");
        }

    }
    if (getLengthBuffer(nts) != getLengthBuffer(quals)) {
        fprintf(stderr,"Sequence and quals differ in length at line %lld: %llu != %llu\n%s\n%s\n",
            (lld) fqr->line, (llu) getLengthBuffer(nts), (llu) getLengthBuffer(quals), 
            getStartBuffer(nts), getStartBuffer(quals));
    }
    int read_len = getLengthBuffer(nts);
    if (read_len > fqr->max_read_len)
        fqr->max_read_len = read_len;
    // reset to the original Buffer
    fqr->buf = lastBuf;
    assert(getStartBuffer(id)[0] == '@');
    return 1;
}

void close_fq(fq_reader_t fqr)
{
    //fprintf(stdout,"close_fq(%s)\n", fqr->name);
#ifndef NO_GZIP
    if (fqr->gz) {
        gzclose_track(fqr->gz);
        fqr->name[0] = '\0';
        fqr->gz = NULL;
    }
#endif
    if (fqr->f) {
        if (fclose_track(fqr->f) != 0) WARN("Could not fclose_track fastq %s! %s\n", fqr->name, strerror(errno));
        fqr->name[0] = '\0';
        fqr->f = NULL;
    }
}


int64_t estimate_fq(char *fname, int sampledReads, int64_t *estimatedReads, int64_t *estimatedBases)
{
    fq_reader_t fqr = create_fq_reader();
    open_fq(fqr, fname, 0);
    int max = sampledReads;
    int i;
    int64_t bases = 0;
    int64_t startPos = fqr->fpos;
    for(i = 0; i < max; i++) {
         if (fqr->fpos >= fqr->end_read) break;
         if (!get_next_line(fqr)) break;
         char *buf = getStartBuffer(fqr->buf);
         if (buf[0] != '@') {
             long long fpos = fqr->fpos;
             fprintf(stderr,"Improper fastq file (%s:%lld) - First line does not start with '@': %s", fname, fpos, buf);
         }
         if (!get_next_line(fqr)) {
             long long fpos = fqr->fpos;
             fprintf(stderr,"Improper fastq file (%s:%lld) - no second line after %s", fname, fpos, buf);
         }
         int seqLen = getLengthBuffer(fqr->buf);
         bases += seqLen - 1;
         if (!get_next_line(fqr)) {
             long long fpos = fqr->fpos;
             fprintf(stderr,"Improper fastq file (%s:%lld) - no third line after: %s", fname, fpos, buf);
         }
         if (buf[0] != '+') {
             long long fpos = fqr->fpos;
             fprintf(stderr,"Improper fastq file (%s:%lld) - third line does not start with '+': %s", fname, fpos, buf);
         }
         if (!get_next_line(fqr)) {
             long long fpos = fqr->fpos;
             fprintf(stderr,"Improper fastq file (%s:%lld) - no fourth line after: %s", fname, fpos, buf);
         }
         if(seqLen != getLengthBuffer(fqr->buf)) {
             long long fpos = fqr->fpos;
             int qualLen = getLengthBuffer(fqr->buf);
             fprintf(stderr,"Improper fastq file (%s:%lld) - seq and qual lengths mismatch: %d vs %d: %s",
                 fname, fpos, seqLen, qualLen, buf);
         }
    }
    int64_t fileSize = fqr->size;
    int64_t bytesRead = fqr->fpos - startPos;
    if (bytesRead > 0) {
        *estimatedReads = (fileSize * i + bytesRead - 1) / bytesRead;
        *estimatedBases = (fileSize * bases + bytesRead - 1) / bytesRead;
    } else {
        fprintf(stdout, "No reads in %s\n", fname);
        *estimatedReads = 0;
        *estimatedBases = 0;
    }
    destroy_fq_reader(fqr);
    return fileSize;
}
