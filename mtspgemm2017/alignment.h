#include <seqan/sequence.h>
#include <seqan/align.h>
#include <seqan/score.h>
#include <seqan/modifier.h>
#include <seqan/seeds.h>
#include "global.h"
#include "alncommon.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <ctype.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>

using namespace seqan;
using namespace std;

alignmentInfo gabaTest(char const *row, char const *col, int rowStart, int colStart, int kmer_len) {

    /* create config */

    //gaba_t *ctx = gaba_init(&ginit);
    gaba_t *ctx = gaba_init(GABA_PARAMS(
        GABA_SCORE_SIMPLE(2, 3, 5, 1),               /* match award, mismatch penalty, gap open penalty (G_i), and gap extension penalty (G_e) */
        gfa : 2,
        gfb : 2,
        xdrop : 100,
        filter_thresh : 0
    ));

    //char const *a = "\x01\x08\x01\x08\x01\x08";       /* 4-bit encoded "ATATAT" */
    //char const *b = "\x01\x08\x01\x02\x01\x08";       /* 4-bit encoded "ATACAT" */
    char const t[64] = { 0 };                           /* tail array */

    gaba_section_s asec = gaba_build_section(0, (uint32_t)strlen(row), row);
    gaba_section_s bsec = gaba_build_section(2, (uint32_t)strlen(col), col);
    gaba_section_s tail = gaba_build_section(4, 64, t);

    /* create thread-local object */
    gaba_dp_t *dp = gaba_dp_init(ctx);                      /* dp[0] holds a 64-cell-wide context */
    // gaba_dp_t *dp_32 = &dp[_dp_ctx_index(32)];           /* dp[1] and dp[2] are narrower ones */
    // gaba_dp_t *dp_16 = &dp[_dp_ctx_index(16)];

    /* init section pointers */
    gaba_section_s const *ap = &asec, *bp = &bsec;
    gaba_fill_s const *f = gaba_dp_fill_root(dp, /* dp -> &dp[_dp_ctx_index(band_width)] makes the band width selectable */
        ap, rowStart,                                   /* a-side (reference side) sequence and start position */
        bp, colStart,                                   /* b-side (query) */
        UINT32_MAX                                      /* max extension length */
    );

    /* until X-drop condition is detected */
    struct gaba_fill_s const *m = f;                    /* track max */
    while((f->status & GABA_TERM) == 0) {
        if(f->status & GABA_UPDATE_A) { ap = &tail; }   /* substitute the pointer by the tail section's if it reached the end */
        if(f->status & GABA_UPDATE_B) { bp = &tail; }

        f = gaba_dp_fill(dp, f, ap, bp, UINT32_MAX);    /* extend the banded matrix */
        m = f->max > m->max ? f : m;                    /* swap if maximum score was updated */
    }

    struct gaba_alignment_s *r = gaba_dp_trace(dp,
        m,                                              /* section with the max */
        NULL                                            /* custom allocator: see struct gaba_alloc_s in gaba.h */
    );

    // Struct gaba_segment_s const *seg is a field of gaba_alignment_s
    alignmentInfo result; // defined in global.h
    result.score = r->score;
    result.apos = r->seg->apos; // row pos
    result.bpos = r->seg->bpos; // col pos
    result.alen = r->seg->alen; // row alignment length
    result.blen = r->seg->blen; // col alignment length

    /* If you need to align the reverse-complement of a sequence you have on memory, 
    you do not need to explicitly construct the reverse-complemented sequence before calling 
    the fill-in functions. Passing a pointer mirrored at GABA_EOU (end-of-userland) will 
    fetch the sequence at the original location in the reverse-complemented manner. 
    gaba_mirror(sequence, strlen(sequence)) will create a pointer to the head of the 
    reverse-complemented sequence (or the tail of the original sequence). */
    /* clean up */
    gaba_dp_res_free(dp, r);
    gaba_dp_clean(dp);
    gaba_clean(ctx);
    return result;
}

/**
 * @brief seqanAlOne does the seed-and-extend alignment
 * when ony one shared k-mer exists
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the k-mer on the first read
 * @param j is the starting position of the k-mer on the second read
 * @param dropFactor
 * @return alignment score and extended seed
 */

seqAnResult seqanAlOne(std::string & row, std::string & col, int rlen, int i, int j, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+kmer_len, j+kmer_len);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        i = rlen-i-kmer_len;
        
        setBeginPositionH(seed1, i);
        setBeginPositionV(seed1, j);
        setEndPositionH(seed1, i+kmer_len);
        setEndPositionV(seed1, j+kmer_len);

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    } 

    longestExtensionScore.score = longestExtensionTemp;
    longestExtensionScore.seed = seed1;
    longestExtensionScore.strand = strand;
    return longestExtensionScore;
}
/**
 * @brief seqanAlGen does the seed-and-extend alignment
 * when two shared k-mers exist
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param i is the starting position of the first k-mer on the first read
 * @param j is the starting position of the first k-mer on the second read
 * @param l is the starting position of the second k-mer on the first read
 * @param m is the starting position of the second k-mer on the second read
 * @param dropFactor
 * @return alignment score and extended seed
 */
seqAnResult seqanAlGen(std::string & row, std::string & col, int rlen, int i, int j, int l, int m, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    std::pair<int,int> longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    TSeed seed1(i, j, i+kmer_len, j+kmer_len);
    TSeed seed2(l, m, l+kmer_len, m+kmer_len);
    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);

        i = rlen-i-kmer_len;
        l = rlen-l-kmer_len;

        setBeginPositionH(seed1, i);
        setBeginPositionV(seed1, j);
        setEndPositionH(seed1, i+kmer_len);
        setEndPositionV(seed1, j+kmer_len);

        setBeginPositionH(seed2, l);
        setBeginPositionV(seed2, m);
        setEndPositionH(seed2, l+kmer_len);
        setEndPositionV(seed2, m+kmer_len);

        /* Perform match extension */
        longestExtensionTemp.first = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp.first = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
        longestExtensionTemp.second = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    }

    longestExtensionScore.strand = strand;

    if(longestExtensionTemp.first > longestExtensionTemp.second)
    {
        longestExtensionScore.score = longestExtensionTemp.first;
        longestExtensionScore.seed = seed1;
    }
    else
    {
        longestExtensionScore.score = longestExtensionTemp.second;
        longestExtensionScore.seed = seed2;
    } 

    return longestExtensionScore;
}
/**
 * @brief seqanAlOneAllKmer does the seed-and-extend alignment
 * when only one shared k-mer exists
 * @param row
 * @param col
 * @param rlen is the length of the row sequence
 * @param vpos is a vector containing the starting position of k-mer
 * @param dropFactor
 * @return alignment score and extended seed
 * DO NOT USE IT UNTIL THE SEED POSITION DEFINITION WILL BE CORRECTED AND UPDATED
 */
seqAnResult seqanAlOneAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int longestExtensionTemp;
    seqAnResult longestExtensionScore;

    seqH = row;
    seqV = col;

    std::vector<std::pair<int,int>>::iterator it;
    it = vpos.begin();

    TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));

    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        it->first = rlen-it->first-kmer_len;

        setBeginPositionH(seed, it->first);
        setBeginPositionV(seed, it->second);
        setEndPositionH(seed, it->first+kmer_len);
        setEndPositionV(seed, it->second+kmer_len);

        //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len); // In this way the seed is not globally updated

        /* Perform match extension */
        longestExtensionTemp = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());

    } else
    {
        strand = 'n';
        longestExtensionTemp = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
    } 

    longestExtensionScore.score = longestExtensionTemp;
    longestExtensionScore.seed = seed;
    longestExtensionScore.strand = strand;

    return longestExtensionScore;
}

/**
 * @brief seqanAlGenAllKmer does the seed-and-extend alignment
 * when shared k-mers > 1
 * @param row
 * @param col
 * @param rlen
 * @param vpos
 * @param dropFactor
 * @return alignment score and extended seed
 * DO NOT USE IT UNTIL THE SEED POSITION DEFINITION WILL BE CORRECTED AND UPDATED
 */
seqAnResult seqanAlGenAllKmer(std::string & row, std::string & col, int rlen, std::vector<std::pair<int,int>> vpos, int dropFactor, int kmer_len) {

    Score<int, Simple> scoringScheme(1, -1, -1);

    Dna5String seqH; 
    Dna5String seqV; 
    Dna5String seedH;
    Dna5String seedV;
    string strand;
    int tempScore;
    seqAnResult longestExtensionScore;
    std::vector<std::pair<int,int>>::iterator fit;

    seqH = row;
    seqV = col;
    longestExtensionScore.score = 0;
    fit = vpos.begin();

    TSeed seed(fit->first, fit->second, fit->first+kmer_len, fit->second+kmer_len);
    seedH = infix(seqH, beginPositionH(seed), endPositionH(seed));
    seedV = infix(seqV, beginPositionV(seed), endPositionV(seed));
    
    Dna5StringReverseComplement twin(seedH);

    if(twin == seedV)
    {   /* Reverse seq */
        strand = 'c';
        Dna5StringReverseComplement twinRead(seqH);
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Update position on reversed seq */
            it->first = rlen-it->first-kmer_len;
            /* Seed unpdate */
            setBeginPositionH(seed, it->first);
            setBeginPositionV(seed, it->second);
            setEndPositionH(seed, it->first+kmer_len);
            setEndPositionV(seed, it->second+kmer_len);
            //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len); // In this way the seed is not globally updated
            /* Perform match extension */
            tempScore = extendSeed(seed, twinRead, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.score)
            {
                longestExtensionScore.score = tempScore;
                longestExtensionScore.seed = seed;
            }
        }
    } 
    else
    {
        strand = 'n';
        for(std::vector<std::pair<int,int>>::iterator it=vpos.begin(); it!=vpos.end(); ++it)
        {   /* Seed update */
            setBeginPositionH(seed, it->first);
            setBeginPositionV(seed, it->second);
            setEndPositionH(seed, it->first+kmer_len);
            setEndPositionV(seed, it->second+kmer_len);
            //TSeed seed(it->first, it->second, it->first+kmer_len, it->second+kmer_len);
            /* Perform match extension */
            tempScore = extendSeed(seed, seqH, seqV, EXTEND_BOTH, scoringScheme, dropFactor, GappedXDrop());
            /* Keep the best score */
            if(tempScore > longestExtensionScore.score)
            {
                longestExtensionScore.score = tempScore;
                longestExtensionScore.seed = seed; 
            }
        }  
    }
    longestExtensionScore.strand = strand;

    return longestExtensionScore;
}
