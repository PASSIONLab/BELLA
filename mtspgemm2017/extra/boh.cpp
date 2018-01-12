#ifdef _SEQAN
        Dna5String seqH; 
        Dna5String seqV; 
        Dna5String seedH;
        Dna5String seedV;
        int64_t longestExtensionScore1;
        int64_t longestExtensionScore2;

        Score<int, Simple> scoringScheme(1, -1, -1);

        #pragma omp parallel for private(myBatch, seedH, seedV, seqH, seqV, longestExtensionScore1, longestExtensionScore2) shared(colStart,colEnd,numCols, reads)
        for(int b = 0; b < numBlocks+1; ++b) 
        { 
            vector<IT> * RowIdsofC = new vector<IT>[numCols[b]];  // row ids for each column of C (bunch of cols)
            vector<FT> * ValuesofC = new vector<FT>[numCols[b]];  // values for each column of C (bunch of cols)
            
            IT * colptr = new IT[numCols[b]+1];
            colptr[0] = 0;
           
            LocalSpGEMM(colStart[b], colEnd[b], numCols[b], A, B, multop, addop, RowIdsofC, ValuesofC);
     
            int k=0;
            for(int i=0; i<numCols[b]; ++i) // for all edge lists (do NOT parallelize)
            {
                colptr[i+1] = colptr[i] + RowIdsofC[k].size();
                ++k;
            }
          
            IT * rowids = new IT[colptr[numCols[b]]];
            FT * values = new FT[colptr[numCols[b]]];
           
            k=0;
            for(int i=0; i<numCols[b]; ++i) // combine step
            {
                copy(RowIdsofC[k].begin(), RowIdsofC[k].end(), rowids + colptr[i]);
                copy(ValuesofC[k].begin(), ValuesofC[k].end(), values + colptr[i]);
                ++k;
            }

            delete [] RowIdsofC;
            delete [] ValuesofC;

            /* Local Alignment before write on stringstream */
            for(int i=0; i<numCols[b]; ++i) 
            {
                for(int j=colptr[i]; j<colptr[i+1]; ++j) 
                {
                    //myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << endl;
                    if(values[j]->count == 1)
                    {          
                        seqH = reads[rowids[j]].seq;
                        seqV = reads[i+colStart[b]].seq;

                        Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                        seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
                        seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

                        Dna5StringReverseComplement twin(seedH);

                        if(twin == seedV)
                        {
                            Dna5StringReverseComplement twinRead(seqH);
                            values[j]->pos[0] = reads[rowids[j]].seq.length()-values[j]->pos[0]-KMER_LENGTH;
                            Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);

                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                
                        } else
                        {
                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
        
                        }

                        if(longestExtensionScore1 > 1200)
                        {
            
                            myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore1 << ' ' << beginPositionV(seed1) << ' ' << 
                                endPositionV(seed1) << ' ' << reads[i+colStart[b]].seq.length() << ' ' << beginPositionH(seed1) << ' ' << endPositionH(seed1) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                        }
                    } 
                    else if(values[j]->count == 2)
                    {       

                        seqH = reads[rowids[j]].seq;
                        seqV = reads[i+colStart[b]].seq;

                        Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                        Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);
                        seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
                        seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

                        Dna5StringReverseComplement twin(seedH);

                        if(twin == seedV)
                        {
                            Dna5StringReverseComplement twinRead(seqH);

                            values[j]->pos[0] = reads[rowids[j]].seq.length()-values[j]->pos[0]-KMER_LENGTH;
                            values[j]->pos[2] = reads[rowids[j]].seq.length()-values[j]->pos[2]-KMER_LENGTH;

                            Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                            Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);

                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                            longestExtensionScore2 = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());

                        } else
                        {
                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                            longestExtensionScore2 = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                        
                        }

                        if(longestExtensionScore1 > 1200 || longestExtensionScore2 > 1200)
                        {
                            if(longestExtensionScore1 > longestExtensionScore2)
                            {
                                myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore1 << ' ' << beginPositionV(seed1) << ' ' << 
                                endPositionV(seed1) << ' ' << reads[i+colStart[b]].seq.length() << ' ' << beginPositionH(seed1) << ' ' << endPositionH(seed1) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                            } else
                            {
                                myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore2 << ' ' << beginPositionV(seed2) << ' ' << 
                                endPositionV(seed2) << ' ' << reads[i+colStart[b]].seq.length() << ' ' << beginPositionH(seed2) << ' ' << endPositionH(seed2) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                            }
                        }
                    } else
                    {
                        seqH = reads[rowids[j]].seq;
                        seqV = reads[i+colStart[b]].seq;

                        Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                        Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);
                        seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
                        seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

                        Dna5StringReverseComplement twin(seedH);

                        if(twin == seedV)
                        {
                            Dna5StringReverseComplement twinRead(seqH);

                            values[j]->pos[0] = reads[rowids[j]].seq.length()-values[j]->pos[0]-KMER_LENGTH;
                            values[j]->pos[2] = reads[rowids[j]].seq.length()-values[j]->pos[2]-KMER_LENGTH;

                            Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                            Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);

                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                            longestExtensionScore2 = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());

                        } else
                        {
                            /* Perform match extension */
                            longestExtensionScore1 = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                            longestExtensionScore2 = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, 5, GappedXDrop());
                        
                        }

                        if(longestExtensionScore1 > 1200 || longestExtensionScore2 > 1200)
                        {
                            if(longestExtensionScore1 > longestExtensionScore2)
                            {
                                myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore1 << ' ' << beginPositionV(seed1) << ' ' << 
                                endPositionV(seed1) << ' ' << reads[i+colStart[b]].seq.length() << ' ' << beginPositionH(seed1) << ' ' << endPositionH(seed1) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                            } else
                            {
                                myBatch << reads[i+colStart[b]].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore2 << ' ' << beginPositionV(seed2) << ' ' << 
                                endPositionV(seed2) << ' ' << reads[i+colStart[b]].seq.length() << ' ' << beginPositionH(seed2) << ' ' << endPositionH(seed2) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                            }
                        }
                    }
                }
            }

            delete [] colptr;
            delete [] rowids;
            delete [] values;

            #pragma omp critical
            {
                writeToFile(myBatch, "ovls.bella");
                myBatch.str(std::string());
            }
        }
        delete [] colStart;
        delete [] colEnd;
        delete [] numCols;
    } else
    {
        int colStart = 0;
        int colEnd = B.cols;
        int numCols = B.cols;
        vector<IT> * RowIdsofC = new vector<IT>[B.cols];      // row ids for each column of C
        vector<FT> * ValuesofC = new vector<FT>[B.cols];      // values for each column of C

        LocalSpGEMM(colStart, colEnd, numCols, A, B, multop, addop, RowIdsofC, ValuesofC);

        IT * colptr = new IT[B.cols+1];
        colptr[0] = 0;
        std::stringstream myBatch;

        for(int i=0; i < B.cols; ++i)        // for all edge lists (do NOT parallelize)
            colptr[i+1] = colptr[i] + RowIdsofC[i].size();

        IT * rowids = new IT[colptr[B.cols]];
        FT * values = new FT[colptr[B.cols]];
        
        for(int i=0; i < B.cols; ++i)        // for all edge lists (do NOT parallelize)
        {
            copy(RowIdsofC[i].begin(), RowIdsofC[i].end(), rowids + colptr[i]);
            copy(ValuesofC[i].begin(), ValuesofC[i].end(), values + colptr[i]);
        }

        delete [] RowIdsofC;
        delete [] ValuesofC;

        Dna5String seqH; 
        Dna5String seqV; 
        Dna5String seedH;
        Dna5String seedV;
        int64_t longestExtensionScore1;
        int64_t longestExtensionScore2;

        std::stringstream myBatchSingle;
        Score<int, Simple> scoringScheme(1, -1, -1);

        for(int i=0; i<B.cols; ++i) 
        {
            for(int j=colptr[i]; j<colptr[i+1]; ++j) 
            {
                if(values[j]->count == 1)
                {          
                    seqH = reads[rowids[j]].seq;
                    seqV = reads[i].seq;

                    Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
                    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

                    Dna5StringReverseComplement twin(seedH);

                    if(twin == seedV)
                    {
                        Dna5StringReverseComplement twinRead(seqH);
                        values[j]->pos[0] = reads[rowids[j]].seq.length()-values[j]->pos[0]-KMER_LENGTH;
                        Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);

                        /* Perform match extension */
                        longestExtensionScore1 = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
                
                    } else
                    {
                        /* Perform match extension */
                        longestExtensionScore1 = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
        
                    }

                    if(longestExtensionScore1 > 800)
                    {
                            myBatchSingle << reads[i].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore1 << ' ' << beginPositionV(seed1) << ' ' << 
                                endPositionV(seed1) << ' ' << reads[i].seq.length() << ' ' << beginPositionH(seed1) << ' ' << endPositionH(seed1) <<
                                    ' ' << reads[rowids[j]].seq.length() << endl;
                    }
                } 
                else if(values[j]->count > 1)
                {       

                    seqH = reads[rowids[j]].seq;
                    seqV = reads[i].seq;

                    Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                    Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);
                    seedH = infix(seqH, beginPositionH(seed1), endPositionH(seed1));
                    seedV = infix(seqV, beginPositionV(seed1), endPositionV(seed1));

                    Dna5StringReverseComplement twin(seedH);

                    if(twin == seedV)
                    {
                        Dna5StringReverseComplement twinRead(seqH);

                        values[j]->pos[0] = reads[rowids[j]].seq.length()-values[j]->pos[0]-KMER_LENGTH;
                        values[j]->pos[2] = reads[rowids[j]].seq.length()-values[j]->pos[2]-KMER_LENGTH;

                        Seed<Simple> seed1(values[j]->pos[0], values[j]->pos[1], values[j]->pos[0]+KMER_LENGTH, values[j]->pos[1]+KMER_LENGTH);
                        Seed<Simple> seed2(values[j]->pos[2], values[j]->pos[3], values[j]->pos[2]+KMER_LENGTH, values[j]->pos[3]+KMER_LENGTH);

                        /* Perform match extension */
                        longestExtensionScore1 = extendSeed(seed1, twinRead, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
                        longestExtensionScore2 = extendSeed(seed2, twinRead, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());

                    } else
                    {
                        /* Perform match extension */
                        longestExtensionScore1 = extendSeed(seed1, seqH, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
                        longestExtensionScore2 = extendSeed(seed2, seqH, seqV, EXTEND_BOTH, scoringScheme, 3, GappedXDrop());
                    
                    }

                    if(longestExtensionScore1 > 800 || longestExtensionScore2 > 800)
                    {
                        if(longestExtensionScore1 > longestExtensionScore2)
                        {
                            myBatchSingle << reads[i].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore1 << ' ' << beginPositionV(seed1) << ' ' << 
                            endPositionV(seed1) << ' ' << reads[i].seq.length() << ' ' << beginPositionH(seed1) << ' ' << endPositionH(seed1) <<
                                ' ' << reads[rowids[j]].seq.length() << endl;
                        } else
                        {
                            myBatchSingle << reads[i].nametag << ' ' << reads[rowids[j]].nametag << ' ' << values[j]->count << ' ' << longestExtensionScore2 << ' ' << beginPositionV(seed2) << ' ' << 
                            endPositionV(seed2) << ' ' << reads[i].seq.length() << ' ' << beginPositionH(seed2) << ' ' << endPositionH(seed2) <<
                                ' ' << reads[rowids[j]].seq.length() << endl;
                        }
                    }
                }
            }
        }
        
        delete [] rowids;
        delete [] values;
        delete [] colptr;

        writeToFile(myBatchSingle, "ovls.bella");
        myBatchSingle.str(std::string());
    }
    #endif