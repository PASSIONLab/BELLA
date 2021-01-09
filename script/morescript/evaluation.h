/**
* Evaluation code to compare long read-to-read aligner recall/precision
* 
* @author: Giulia Guidi
*
*/
#include "IntervalTree.h"
#include <omp.h>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
#include <sys/types.h> 
#include <sys/stat.h> 
#include <math.h>
#include <limits.h>
#include <bitset>
#include <unordered_map>
#include <map>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <ctype.h> 
#include <sstream>
#include <set>
#include <memory>
#include <typeinfo>

using namespace std;
#define MOREINFO

struct refInfo {

    string ref;
    string read;
    int start;
    int end;
};

struct readInfo {

    string ref;
    int start;
    int end;
};

#ifdef MOREINFO
struct myStruct {

    int state;
    string aStart;
    string aEnd;
    string aLen;
    string bStart;
    string bEnd;
    string bLen;
    string skmer;
    string score;
};
#endif

typedef vector<refInfo> vectInfo;

int estimeOvl(int colStart, int colEnd, int colLen, int rowStart, int rowEnd, int rowLen)
{
    int diffCol = colEnd-colStart;
    int diffRow = rowEnd-rowStart;
    int minLeft = min(colStart, rowStart);
    int minRight = min(colLen-colEnd, rowLen-rowEnd);

    int result = minLeft+minRight+(diffCol+diffRow)/2;
    return result;
}

/**
 * @brief trueOv omputes the number of true overlapping reads
 * @param truthInfo
 * @param simulated (default false)
 * @param minOvl
 * @return number of true overlapping pairs
 */
double trueOv(map<string,vectInfo> & truthInfo, bool simulated, int minOvl) 
{
    vector<Interval<std::string>> intervals;
    vector<Interval<std::string>> queries;
    vector<Interval<std::string>>::iterator q;
    map<string,vectInfo>::iterator key; // outer iterator
    vectInfo::iterator it; // inner iterator
    double trueOvls = 0;

    for(key = truthInfo.begin(); key != truthInfo.end(); ++key)
    {
         for(it = key->second.begin(); it != key->second.end(); ++it)
         {
             intervals.push_back(Interval<string>(it->start, it->end, it->read));
             queries.push_back(Interval<string>(it->start, it->end, it->read));
         }

         IntervalTree<string> tree;
         vector<size_t> treeCount;
         
         tree = IntervalTree<string>(intervals); // reference 
         
         for (q = queries.begin(); q != queries.end(); ++q) // tree search for a given reference
         {
             vector<Interval<string>> results;
             tree.findOverlapping(q->start, q->stop, q->value, results, minOvl);
             treeCount.push_back(results.size());
         }
         
         for(size_t t = 0; t < treeCount.size(); ++t) // count for a given reference
         { 
             trueOvls = trueOvls + (double)treeCount[t];  // cumulative
         }

         intervals.clear();
         queries.clear();
    }  
    return trueOvls;
}

/**
 * @brief computeLength computes the overlap length between
 * potential overlapping reads pairs
 * @param readMap
 * @param colName
 * @param rowName
 * @return alignment length
 */

typedef multimap<string,readInfo>::iterator MMAPIterator;
int computeLength(multimap<string,readInfo> & readMap, string & colName, string & rowName) 
{
    int length = 0, maxlength = 0;

    MMAPIterator checkCol = readMap.find(colName); // col name
    MMAPIterator checkRow = readMap.find(rowName); // row name 

    if(checkCol != readMap.end() && checkRow != readMap.end()) // needed as handling real dataset the aligned reads in sam file could be != the original number of reads
    {
        // It returns a pair representing the range of elements with key equal to the one we want
        pair<MMAPIterator, MMAPIterator> colRange = readMap.equal_range(colName);
        pair<MMAPIterator, MMAPIterator> rowRange = readMap.equal_range(rowName);
        // Iterate over the range
        for (MMAPIterator itC = colRange.first; itC != colRange.second; itC++)
        {
            for (MMAPIterator itR = rowRange.first; itR != rowRange.second; itR++)
            {
                if(itC->second.ref == itR->second.ref)
                {   
                    if(itC->second.start < itR->second.start) 
                    {
                        if(itC->second.end > itR->second.start)
                            length = min((itC->second.end - itR->second.start), (itR->second.end - itR->second.start));
                    }
                else if(itC->second.start > itR->second.start) 
                {
                    if(itR->second.end > itC->second.start)
                        length = min((itR->second.end - itC->second.start), (itC->second.end - itC->second.start));
                } 
                else length = min((itR->second.end - itC->second.start), (itC->second.end - itC->second.start)); 
                }

                if(length > maxlength)
                    maxlength = length; 
            }
        }
    }
    return maxlength;
}

/**
 * @brief benchmarkingAl retrives recall/precision values
 * @param groundtruth (input file)
 * @param bella (file containing bella alignments)
 * @param minimap (file containing minimap overlaps)
 * @param mhap (file containing mhap overlaps)
 * @param blasr (file containing blasr alignments)
 * @param daligner (file containing daligner alignments)
 * @param simulated (default false)
 * @param minOvl
 */
void benchmarkingAl(ifstream & groundtruth, ifstream & bella, ifstream & minimap, ifstream & mhap, 
    ifstream & blasr, ifstream & daligner, bool simulated, int minOvl) 
{
    map<string,vectInfo> isInThere;
    map<string,vectInfo>::iterator iter;
    multimap<string,readInfo> readMap;

#ifdef MOREINFO
    map<pair<string,string>, myStruct> checkBella;
    map<pair<string,string>, myStruct>::iterator extraIt;
    myStruct extraInfo;
#else
    map<pair<string,string>, int> checkBella;
#endif
    map<pair<string,string>, int> checkMinimap;
    map<pair<string,string>, int> checkMhap;
    map<pair<string,string>, int> checkBlasr;
    map<pair<string,string>, int> checkDaligner;
    map<pair<string,string>, int>::iterator it;
    
    int alignment_length;
    
    double ovlsbella = 0, truebella = 0;
    double ovlsminimap = 0, trueminimap = 0;
    double ovlsmhap = 0, truemhap = 0;
    double ovlsblasr = 0, trueblasr = 0;
    double ovlsdal = 0, truedal = 0;

    cout << "\nbuilding the ground truth" << endl;

    if(simulated)
    {
        if(groundtruth.is_open())
        {
            string ref;
            string prev;
            string read;
            string dontcare1;
            string dontcare2;
            refInfo ovlInfo;
            readInfo perRead;
            int start;
            int end;

            while(groundtruth >> ref >> start >> end >> read >> dontcare1 >> dontcare2)
            {
                perRead.ref = ref;
                perRead.start = start;
                perRead.end = end;
                readMap.insert(make_pair("@"+read,perRead));

                ovlInfo.ref = ref;
                ovlInfo.read = "@" + read;
                ovlInfo.start = start;
                ovlInfo.end = end;
                
                iter = isInThere.find(ref);
                if(iter == isInThere.end())
                {
                    vectInfo temp;
                    temp.push_back(ovlInfo); // all the element of a chromosome
                    isInThere.insert(map<string,vectInfo>::value_type(ref,temp));
                }
                else
                {
                    iter->second.push_back(ovlInfo);
                    isInThere[ref] = iter->second;
                }
            }
            //isInThere.push_back(vectOf); // insert the last chromosome
            cout << "num reads: " << readMap.size() << endl;
            cout << "num sub-references: " << isInThere.size() << endl;
            //for(int i = 0; i < truthInfo.size(); ++i)
            //    cout << "size chromosome: " << truthInfo.at(i).size() << endl;

        } else cout << "Error opening the ground truth file" << endl;
    }
    else // from sam file
    {
        if(groundtruth.is_open())
        {
            refInfo ovlInfo;
            readInfo perRead;
            string prev;
            string read;
            string ref;
            int start;
            int end;
    
            while(groundtruth >> ref >> read >> start >> end)
            {
                perRead.ref = ref;
                perRead.start = start;
                perRead.end = end;
                readMap.insert(make_pair("@"+read,perRead));

                ovlInfo.ref = ref;
                ovlInfo.read = "@" + read;
                ovlInfo.start = start;
                ovlInfo.end = end;
                
                iter = isInThere.find(ref);
                if(iter == isInThere.end())
                {
                    vectInfo temp;
                    temp.push_back(ovlInfo); // all the element of a chromosome
                    isInThere.insert(map<string,vectInfo>::value_type(ref,temp));
                }
                else
                {
                    iter->second.push_back(ovlInfo);
                    isInThere[ref] = iter->second;
                }
            }
            //isInThere.push_back(vectOf); // insert the last chromosome
            cout << "num reads: " << readMap.size() << endl;
            cout << "num sub-references: " << isInThere.size() << endl;
            //for(int i = 0; i < truthInfo.size(); ++i)
                //cout << "size chromosome: " << truthInfo.at(i).size() << endl;

        } else cout << "Error opening the ground truth file" << endl;
    }
    
    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);

    cout << "computing BELLA recall/precision" << endl;
    if(bella.is_open())
    {
        string line;
        while(getline(bella, line))
        {   
            stringstream lineStream(line);
            string colName, rowName, nkmer, score;
            string colStart, colEnd, colLen, rowStart, rowEnd, rowLen, strand;
            
            getline(lineStream, colName, '\t');
            getline(lineStream, rowName, '\t');
            getline(lineStream, nkmer, '\t');
            getline(lineStream, score, '\t');
            getline(lineStream, strand, '\t');
            getline(lineStream, colStart, '\t');
            getline(lineStream, colEnd, '\t');
            getline(lineStream, colLen, '\t');
            getline(lineStream, rowStart, '\t');
            getline(lineStream, rowEnd, '\t');
            getline(lineStream, rowLen, '\t');

            if(colName != rowName) // to be sure to not count self aligned pairs
            {
                // value == 0 ==> pass alignment threshold
                // value == 1 ==> doesn't pass alignment threshold AND pass overlap estimate
                // value == 2 ==> doesn't pass alignment threshold AND doesn't pass overlap estimate
            #ifdef MOREINFO
                extraIt = checkBella.find(make_pair(colName, rowName));
            #else
                it = checkBella.find(make_pair(colName, rowName));
            #endif
            #ifdef MOREINFO
                if(extraIt == checkBella.end())
            #else
                if(it == checkBella.end())
            #endif
                {
                    ovlsbella++;
                    // Compute the overlap length between potential overlapping reads pairs 
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                #ifdef MOREINFO
                        extraInfo.aStart = colStart;
                        extraInfo.aEnd = colEnd;
                        extraInfo.aLen = colLen;
                        extraInfo.bStart = rowStart;
                        extraInfo.bEnd = rowEnd;
                        extraInfo.bLen = rowLen;
                        extraInfo.skmer = nkmer;
                        extraInfo.score = score;
                        extraInfo.state = 0;
                        checkBella.insert(make_pair(make_pair(colName, rowName), extraInfo));
                #else
                        checkBella.insert(make_pair(make_pair(colName, rowName), 0));
                #endif
                        truebella++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                #ifdef MOREINFO
                            extraInfo.aStart = colStart;
                            extraInfo.aEnd = colEnd;
                            extraInfo.aLen = colLen;
                            extraInfo.bStart = rowStart;
                            extraInfo.bEnd = rowEnd;
                            extraInfo.bLen = rowLen;
                            extraInfo.skmer = nkmer;
                            extraInfo.score = score;
                            extraInfo.state = 2;
                            checkBella.insert(make_pair(make_pair(colName, rowName), extraInfo));
                #else
                            checkBella.insert(make_pair(make_pair(colName, rowName), 2));
                #endif
                            ovlsbella--;
                        } else
                        {
                #ifdef MOREINFO
                            extraInfo.aStart = colStart;
                            extraInfo.aEnd = colEnd;
                            extraInfo.aLen = colLen;
                            extraInfo.bStart = rowStart;
                            extraInfo.bEnd = rowEnd;
                            extraInfo.bLen = rowLen;
                            extraInfo.skmer = nkmer;
                            extraInfo.score = score;
                            extraInfo.state = 1;
                            checkBella.insert(make_pair(make_pair(colName, rowName), extraInfo));
                #else
                            checkBella.insert(make_pair(make_pair(colName, rowName), 1));
                #endif
                        } 
                    }
                }
            #ifdef MOREINFO
                else if(extraIt->second.state == 1)
            #else
                else if(it->second == 1) 
            #endif
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                #ifdef MOREINFO
                        checkBella[make_pair(colName, rowName)].state = 0;
                #else
                        checkBella[make_pair(colName, rowName)] = 0;
                #endif
                        truebella++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                #ifdef MOREINFO
                            checkBella[make_pair(colName, rowName)].state = 2;
                #else
                            checkBella[make_pair(colName, rowName)] = 2;
                #endif
                            ovlsbella--;
                        } else
                        {
                #ifdef MOREINFO
                            checkBella[make_pair(colName, rowName)].state = 1;
                #else
                            checkBella[make_pair(colName, rowName)] = 1;
                #endif
                        } 
                    }
                }
            #ifdef MOREINFO
                else if(extraIt->second.state == 2)
            #else
                else if(it->second == 2) 
            #endif
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                #ifdef MOREINFO
                        checkBella[make_pair(colName, rowName)].state = 0;
                #else
                        checkBella[make_pair(colName, rowName)] = 0;
                #endif
                        truebella++;
                        ovlsbella++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                #ifdef MOREINFO
                            checkBella[make_pair(colName, rowName)].state = 2;
                #else
                            checkBella[make_pair(colName, rowName)] = 2;
                #endif
                        }
                        else
                        {
                #ifdef MOREINFO
                            checkBella[make_pair(colName, rowName)].state = 1;
                #else            
                            checkBella[make_pair(colName, rowName)] = 1;
                #endif
                            ovlsbella++;
                        } 
                    }
                }
            }
        }
    } 

#ifdef MOREINFO
    ofstream output("postEvaluation_ecoli30x_full.out");
    if(output.is_open())
    {
        // State, shared k-mer, score, nameA, nameB, startA, endA, lenA, startB, endB, lenB
        for(extraIt = checkBella.begin(); extraIt != checkBella.end(); ++extraIt)
        {
            output << extraIt->second.state << '\t' << extraIt->second.skmer << '\t' << extraIt->second.score << '\t' <<  extraIt->first.first << '\t' << extraIt->first.second << '\t' << extraIt->second.aStart 
                << '\t' << extraIt->second.aEnd << '\t' << extraIt->second.aLen << '\t' << extraIt->second.bStart << '\t' << extraIt->second.bEnd 
                    << '\t' << extraIt->second.bLen << endl; 
        }
    }
    output.close();
#endif
    cout << "computing Minimap recall/precision" << endl;
    if(minimap.is_open())
    {
        string line;
        while(getline(minimap, line))
        {
            stringstream lineStream(line);
            string colName, rowName, colStart, colEnd, rowLen, colLen;
            string strand, rowStart, rowEnd;
    
            getline(lineStream, colName, '\t' );
            getline(lineStream, colLen, '\t' );
            getline(lineStream, colStart, '\t' );
            getline(lineStream, colEnd, '\t' );
            getline(lineStream, strand, '\t' );
            getline(lineStream, rowName, '\t' );
            getline(lineStream, rowLen, '\t' );
            getline(lineStream, rowStart, '\t' );
            getline(lineStream, rowEnd, '\t' );
    
            colName = "@" + colName;
            rowName = "@" + rowName;

            if(colName != rowName) // to be sure to not count self aligned pairs
            {
                // value == 0 ==> pass alignment threshold
                // value == 1 ==> doesn't pass alignment threshold AND pass overlap estimate
                // value == 2 ==> doesn't pass alignment threshold AND doesn't pass overlap estimate
                it = checkMinimap.find(make_pair(colName, rowName));
                if(it == checkMinimap.end())
                {
                    ovlsminimap++;
                    // Compute the overlap length between potential overlapping reads pairs 
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMinimap.insert(make_pair(make_pair(colName, rowName), 0));
                        trueminimap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkMinimap.insert(make_pair(make_pair(colName, rowName), 2));
                            ovlsminimap--;
                        } else checkMinimap.insert(make_pair(make_pair(colName, rowName), 1));
                    }
                }
                else if(it->second == 1)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMinimap[make_pair(colName, rowName)] = 0;
                        trueminimap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkMinimap[make_pair(colName, rowName)] = 2;
                            ovlsminimap--;
                        } else checkMinimap[make_pair(colName, rowName)] = 1;
                    }
                }
                else if(it->second == 2)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMinimap[make_pair(colName, rowName)] = 0;
                        trueminimap++;
                        ovlsminimap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                            checkMinimap[make_pair(colName, rowName)] = 2;
                        else
                        {
                            checkMinimap[make_pair(colName, rowName)] = 1;
                            ovlsminimap++;
                        } 
                    }
                }
            }
        }
    }

    cout << "computing MHAP sensitive recall/precision" << endl;
    if(mhap.is_open())
    {
        string line;
        while(getline(mhap, line))
        {
            stringstream lineStream(line);
            string colName, rowName, error, shared, colStart, colEnd, colLen;
            string rowLen, rowEnd, rowStart, rowStrand, colStrand;

            // [A ID] [B ID] [% error] [# shared min-mers] [0=A fwd, 1=A rc] [A start] [A end] [A length] [0=B fwd, 1=B rc] [B start] [B end] [B length]
            getline(lineStream, colName, ' ');
            getline(lineStream, rowName, ' ');
            getline(lineStream, error, ' ');
            getline(lineStream, shared, ' ');
            getline(lineStream, colStrand, ' ');
            getline(lineStream, colStart, ' ');
            getline(lineStream, colEnd, ' ');
            getline(lineStream, colLen, ' ');
            getline(lineStream, rowStrand, ' ');
            getline(lineStream, rowStart, ' ');
            getline(lineStream, rowEnd, ' ');
            getline(lineStream, rowLen, ' ');

            colName = "@" + colName;
            rowName = "@" + rowName;

            if(colName != rowName) // to be sure to not count self aligned pairs
            {
                // value == 0 ==> pass alignment threshold
                // value == 1 ==> doesn't pass alignment threshold AND pass overlap estimate
                // value == 2 ==> doesn't pass alignment threshold AND doesn't pass overlap estimate
                it = checkMhap.find(make_pair(colName, rowName));
                if(it == checkMhap.end())
                {
                    ovlsmhap++;
                    // Compute the overlap length between potential overlapping reads pairs 
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMhap.insert(make_pair(make_pair(colName, rowName), 0));
                        truemhap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkMhap.insert(make_pair(make_pair(colName, rowName), 2));
                            ovlsmhap--;
                        } else checkMhap.insert(make_pair(make_pair(colName, rowName), 1));
                    }
                }
                else if(it->second == 1)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMhap[make_pair(colName, rowName)] = 0;
                        truemhap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkMhap[make_pair(colName, rowName)] = 2;
                            ovlsmhap--;
                        } else checkMhap[make_pair(colName, rowName)] = 1;
                    }
                }
                else if(it->second == 2)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkMhap[make_pair(colName, rowName)] = 0;
                        truemhap++;
                        ovlsmhap++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                            checkMhap[make_pair(colName, rowName)] = 2;
                        else
                        {
                            checkMhap[make_pair(colName, rowName)] = 1;
                            ovlsmhap++;
                        } 
                    }
                }
            }
        }
    }

    cout << "computing BLASR recall/precision" << endl;
    if(blasr.is_open())
    {
        string line;
        while(getline(blasr, line))
        {
            stringstream lineStream(line);
            string colName, rowName, score, similarity, colStart, colEnd, colLen;
            string rowLen, rowEnd, rowStart, rowStrand, colStrand, mapQV;

            // qName tName score percentSimilarity qStrand qStart qEnd qLength tStrand tStart tEnd tLength mapQV
            getline(lineStream, colName, ' ');
            getline(lineStream, rowName, ' ');
            getline(lineStream, score, ' ');
            getline(lineStream, similarity, ' ');
            getline(lineStream, colStrand, ' ');
            getline(lineStream, colStart, ' ');
            getline(lineStream, colEnd, ' ');
            getline(lineStream, colLen, ' ');
            getline(lineStream, rowStrand, ' ');
            getline(lineStream, rowStart, ' ');
            getline(lineStream, rowEnd, ' ');
            getline(lineStream, rowLen, ' ');
            getline(lineStream, mapQV, ' ');

            colName = "@" + colName;
            rowName = "@" + rowName;

            if(colName != rowName) // to be sure to not count self aligned pairs
            {
                // value == 0 ==> pass alignment threshold
                // value == 1 ==> doesn't pass alignment threshold AND pass overlap estimate
                // value == 2 ==> doesn't pass alignment threshold AND doesn't pass overlap estimate
                it = checkBlasr.find(make_pair(colName, rowName));
                if(it == checkBlasr.end())
                {
                    ovlsblasr++;
                    // Compute the overlap length between potential overlapping reads pairs 
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkBlasr.insert(make_pair(make_pair(colName, rowName), 0));
                        trueblasr++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkBlasr.insert(make_pair(make_pair(colName, rowName), 2));
                            ovlsblasr--;
                        } else checkBlasr.insert(make_pair(make_pair(colName, rowName), 1));
                    }
                }
                else if(it->second == 1)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkBlasr[make_pair(colName, rowName)] = 0;
                        trueblasr++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                        {
                            checkBlasr[make_pair(colName, rowName)] = 2;
                            ovlsblasr--;
                        } else checkBlasr[make_pair(colName, rowName)] = 1;
                    }
                }
                else if(it->second == 2)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkBlasr[make_pair(colName, rowName)] = 0;
                        trueblasr++;
                        ovlsblasr++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), stoi(rowStart), stoi(rowEnd), stoi(rowLen));
                        if (ovlEstime < minOvl)
                            checkBlasr[make_pair(colName, rowName)] = 2;
                        else
                        {
                            checkBlasr[make_pair(colName, rowName)] = 1;
                            ovlsblasr++;
                        } 
                    }
                }
            }
        }
    }

    cout << "computing DALIGNER recall/precision" << endl;
    if(daligner.is_open())
    {
        string line;
        while(getline(daligner, line))
        {   
            //
            // 1  1,865 c   [18,479..20,216] x [ 1,707..0>  (24,451 x 7,283 bps, 19 trace pts)
            //
            // "The display of an LA always begins with a line giving the A-read, then the B-read, 
            // then an indication of orientation (i.e. 'n' for same strand, and 'c' for the opposite strand) 
            // followed by the A-interval and B-interval that are aligned and in parentheses the lengths of 
            // the two reads and the number of tracepoints in the alignment between them. In particular, note 
            // carefully that when the B-read is in the complement orientation (c), then the B-interval gives the 
            // higher coordinate first, the idea being that one will align from the highest base down to the lowest 
            // base in the descending direction on B, complement the characters as you go. Further note that in the 
            // alignment display the coordinates at the start of each line follow this orientation convention and 
            // give the coordinate of the "tick mark" just left of the first character in each line. It is useful 
            // to know if an interval reaches the end of read, and to signal this we use an angle-bracket <> instead 
            // of a square bracket [], e.g. in the example the B-segment starts at the beginning of the read."
            // 
            // After our parsing Daligner out is something similar to the other tools output format:
            // colName rowName strand S1 E1 S2/E2 S2/E2 L1 L2
            // S2/E2 depends on the strand, if strand == n, will be S2 and then E2, otherwise it will be E2, and then S2
            //

            stringstream lineStream(line);
            string colName, rowName, strand, colStart, colEnd, colLen, rowStart, rowEnd, rowLen;
            int irowStart, irowEnd, irowLen;

            getline(lineStream, colName, ' ');
            getline(lineStream, rowName, ' ');
            getline(lineStream, strand, ' ');
            getline(lineStream, colStart, ' ');
            getline(lineStream, colEnd, ' ');
            getline(lineStream, colLen, ' ');
            getline(lineStream, rowStart, ' '); // check is stoi can parse ',' 
            getline(lineStream, rowEnd, ' ');
            getline(lineStream, rowLen, ' ');
            irowStart = stoi(rowStart);
            irowEnd = stoi(rowEnd);
            irowLen = stoi(rowLen);
            if(strand.compare("c") == 0)
            {
                irowStart = irowLen-irowEnd-1;
                irowEnd = irowLen-irowStart-1;
            }

            colName = "@" + colName;
            rowName = "@" + rowName;

            if(colName != rowName) // to be sure to not count self aligned pairs
            {    
                // value == 0 ==> pass alignment threshold
                // value == 1 ==> doesn't pass alignment threshold AND pass overlap estimate
                // value == 2 ==> doesn't pass alignment threshold AND doesn't pass overlap estimate
                it = checkDaligner.find(make_pair(colName, rowName));
                if(it == checkDaligner.end())
                {
                    ovlsdal++;
                    // Compute the overlap length between potential overlapping reads pairs 
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkDaligner.insert(make_pair(make_pair(colName, rowName), 0));
                        truedal++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), irowStart, irowEnd, irowLen);
                        if (ovlEstime < minOvl)
                        {
                            checkDaligner.insert(make_pair(make_pair(colName, rowName), 2));
                            ovlsdal--;
                        } else checkDaligner.insert(make_pair(make_pair(colName, rowName), 1));
                    }
                }
                else if(it->second == 1)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkDaligner[make_pair(colName, rowName)] = 0;
                        truedal++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), irowStart, irowEnd, irowLen);
                        if (ovlEstime < minOvl)
                        {
                            checkDaligner[make_pair(colName, rowName)] = 2;
                            ovlsdal--;
                        } else checkDaligner[make_pair(colName, rowName)] = 1;
                    }
                }
                else if(it->second == 2)
                {
                    alignment_length = computeLength(readMap, colName, rowName);
                    if(alignment_length >= minOvl)
                    {
                        checkDaligner[make_pair(colName, rowName)] = 0;
                        truedal++;
                        ovlsdal++;
                    }
                    else 
                    {
                        int ovlEstime = estimeOvl(stoi(colStart), stoi(colEnd), stoi(colLen), irowStart, irowEnd, irowLen);
                        if (ovlEstime < minOvl)
                            checkDaligner[make_pair(colName, rowName)] = 2;
                        else
                        {
                            checkDaligner[make_pair(colName, rowName)] = 1;
                            ovlsdal++;
                        } 
                    }
                }
            }
        }
    }

    groundtruth.clear();
    groundtruth.seekg(0, ios::beg);
    double truetruth = trueOv(isInThere, simulated, minOvl);

    bella.close();
    minimap.close();
    mhap.close();
    blasr.close();
    daligner.close();
    groundtruth.close();

    // Ground Truth 
    cout << "true overlapping from ground truth = " << truetruth << "\n" << endl;
    // BELLA Recall and precision 
    cout << "overlapping from BELLA = " << ovlsbella << endl;
    cout << "true overlapping from BELLA = " << truebella << endl;
    cout << "recall BELLA = " << truebella/truetruth << endl;
    cout << "precision BELLA = " << truebella/ovlsbella << "\n" << endl;
    // Minimap Recall and precision 
    // as -S option count overlaps only once 
    // (A ov B, but not B ov A), while all other 
    // (included the ground truth) count all ovls and the self-ovls
    cout << "overlapping from minimap = " << ovlsminimap << endl;
    cout << "true overlapping from minimap = " << trueminimap << endl;
    cout << "recall minimap = " << (trueminimap*2)/truetruth << endl;
    cout << "precision minimap = " << trueminimap/ovlsminimap << "\n" << endl;
    // MHAP Recall and precision 
    cout << "overlapping from MHAP = " << ovlsmhap << endl;
    cout << "true overlapping from MHAP= " << truemhap << endl;
    cout << "recall MHAP = " << truemhap/truetruth << endl;
    cout << "precision MHAP = " << truemhap/ovlsmhap << "\n" << endl;
    // BLASR Recall and precision 
    cout << "overlapping from BLASR = " << ovlsblasr << endl;
    cout << "true overlapping from BLASR = " << trueblasr << endl;
    cout << "recall BLASR = " << trueblasr/truetruth << endl;
    cout << "precision BLASR = " << trueblasr/ovlsblasr << "\n" << endl;
    // DALIGNER Recall and precision 
    cout << "overlapping from DALIGNER = " << ovlsdal << endl;
    cout << "true overlapping from DALIGNER = " << truedal << endl;
    cout << "recall DALIGNER = " << truedal/truetruth << endl;
    cout << "precision DALIGNER = " << truedal/ovlsdal << "\n" << endl;
}
