#ifndef _CSV_H_
#define _CSV_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cassert>

using namespace std;

template<class IT, class FT>
class CSVwriter
{
public:
    CSVwriter(IT & ncols, IT & pcols, vector<IT> * colptr, vector<IT> * rowids, vector<FT> * values, std::ofstream * stream) // Default constructor
    {
    	maxrows = ncols;
    	skiprows = pcols;
    	copy(colptr.begin(), colptr.end(), batchC);
    	copy(rowids.begin(), rowids.end(), batchR);
    	copy(values.begin(), values.end(), batchV);
    	Open(stream);
    	_delimiter = ",";
    }

    bool Close();
	bool Open();
	void WriteBlock();

	~CSVwriter()
	{
		if(colptr.size() > 0)
			delete [] batchC;
		if(rowids.size() > 0)
			delete [] batchR;
		if(values.size() > 0)
			delete [] batchV;
	}

	IT maxrows;
	IT skiprows;
	vector<IT> * batchC;
    vector<IT> * batchR;
    vector<NT> * batchV;
    unsigned char _delimiter; // = ,
}

class CSVreader
{
public:
    CSVreader(IT & ncols, IT & pcols, std::ifstream * stream) // Default constructor
    {
    	maxrows = ncols;
    	skiprows = pcols;
    	batchC = NULL;
    	batchR = NULL;
    	batchV = NULL;
    	Open(stream);
    	_delimiter = ",";
    }

    bool Close();
	bool Open();
	void ReadBlock();

	~CSVreader()
	{
		if(batchC.size() > 0)
			delete [] batchC;
		if(batchR.size() > 0)
			delete [] batchR;
		if(batchV.size() > 0)
			delete [] batchV;
	}

	IT maxrows;
	IT skiprows;
	vector<IT> * batchC;
    vector<IT> * batchR;
    vector<NT> * batchV;
    unsigned char _delimiter; // = ,
}

#include "CSV.cpp"
#endif