#include "CSC.h"
#include "utility.h"
#include "BitMap.h"
#include <algorithm>
#include <numeric>
using namespace std;

// copy constructor
template <class IT, class NT> 
CSC<IT,NT>::CSC (const CSC<IT,NT> & rhs): nnz(rhs.nnz), rows(rhs.rows), cols(rhs.cols)
{
	if(nnz > 0)
	{
		values = new NT[nnz]; 
		rowids = new IT[nnz];   
        copy(rhs.values, rhs.values + nnz, values);
        copy(rhs.rowids, rhs.rowids + nnz, rowids);
	}
	if ( cols > 0)
	{
		colptr	= new IT[cols+1];
        copy(rhs.colptr, rhs.colptr + cols+1, colptr);
	}
}

template <class IT, class NT>
CSC<IT,NT> & CSC<IT,NT>::operator= (const CSC<IT,NT> & rhs) 
{
	if(this != &rhs)		
	{
		if(nnz > 0)	// if the existing object is not empty
		{
			delete [] rowids;   // empty it
			delete [] values;
		}
		if(cols > 0)
		{
			delete [] colptr;
		}

		nnz	= rhs.nnz;
		rows = rhs.rows;
		cols = rhs.cols;
		if(rhs.nnz > 0)	// if the copied object is not empty
		{
			values = new NT[nnz];
			rowids = new IT[nnz];
            copy(rhs.values, rhs.values + nnz, values);
            copy(rhs.rowids, rhs.rowids + nnz, rowids);
		}
		if(rhs.cols > 0)
		{
			colptr	= new IT[cols+1];
            copy(rhs.colptr, rhs.colptr + cols+1, colptr);
		}
	}
	return *this;
}


//! Construct a CSC object from a GTgraph object
//! GTgraph might have parallel edges; this constructor sums them up
//! Assumes a sorted GTgraph (primary key: start)
template <class IT, class NT>
CSC<IT,NT>::CSC(graph & G):nnz(G.m), rows(G.n), cols(G.n)
{
	// graph is like a triples object
	// typedef struct {
        // LONG_T m;
        // LONG_T n;
        // // Arrays of size 'm' storing the edge information
        // // A directed edge 'e' (0 <= e < m) from start[e] to end[e]
        // // had an integer weight w[e] 
        // LONG_T* start;
        // LONG_T* end; 
	// WEIGHT_T* w;
	// } graph; 
	cout << "Graph nnz= " << G.m << " and n=" << G.n << endl;

	vector< Triple<IT,NT> > simpleG;
	vector< pair< pair<IT,IT>,NT> > currCol;
	currCol.push_back(make_pair(make_pair(G.start[0], G.end[0]), G.w[0]));
	for (IT k = 0 ; k < nnz-1 ; ++k)
        {
                if(G.start[k] != G.start[k+1] )
                {
			std::sort(currCol.begin(), currCol.end());
			simpleG.push_back(Triple<IT,NT>(currCol[0].first.first, currCol[0].first.second, currCol[0].second));
			for(int i=0; i< currCol.size()-1; ++i)
			{
				if(currCol[i].first == currCol[i+1].first)
				{
					simpleG.back().val += currCol[i+1].second;
				}
				else
				{	
					simpleG.push_back(Triple<IT,NT>(currCol[i+1].first.first, currCol[i+1].first.second, currCol[i+1].second));
				}
			}
			vector< pair< pair<IT,IT>,NT> >().swap(currCol);
		}
		currCol.push_back(make_pair(make_pair(G.start[k+1], G.end[k+1]), G.w[k+1]));
        }
	// now do the last row
	sort(currCol.begin(), currCol.end());
        simpleG.push_back(Triple<IT,NT>(currCol[0].first.first, currCol[0].first.second, currCol[0].second));
        for(int i=0; i< currCol.size()-1; ++i)
        {
        	if(currCol[i].first == currCol[i+1].first)
               	{
               		simpleG.back().val += currCol[i+1].second;
		}
		else
                {
        	        simpleG.push_back(Triple<IT,NT>(currCol[i+1].first.first, currCol[i+1].first.second, currCol[i+1].second));
                }
        }

	nnz = simpleG.size();
	cout << "[After duplicate merging] Graph nnz= " << nnz << " and n=" << G.n << endl;

	colptr = new IT[cols+1];
    	rowids = new IT[nnz];
    	values = new NT[nnz];

	IT * work = new IT[cols];	// workspace
    	std::fill(work, work+cols, (IT) 0); // initilized to zero
   
    	for (IT k = 0 ; k < nnz ; ++k)
    	{
        	IT tmp =  simpleG[k].col;
        	work [ tmp ]++ ;		// col counts (i.e, w holds the "col difference array")
	}

	if(nnz > 0)
	{
		colptr[cols] = CumulativeSum (work, cols) ;		// cumulative sum of w
        	copy(work, work+cols, colptr);

		IT last;
		for (IT k = 0 ; k < nnz ; ++k)
        	{
			rowids[ last = work[ simpleG[k].col ]++ ]  = simpleG[k].row ;
			values[last] = simpleG[k].val ;
       	 	}
	}
	delete [] work;
}

// Construct a Csc object from an array of "triple"s 
template <class IT, class NT>
CSC<IT,NT>::CSC(Triple<IT,NT> * triples, IT mynnz, IT m, IT n):nnz(mynnz),rows(m),cols(n)
{
	colptr = new IT[cols+1]();
    	rowids = new IT[nnz];
    	values = new NT[nnz];
	vector< pair<IT,NT> > tosort (nnz);

	IT * work = new IT[cols];	// workspace
    std::fill(work, work+cols, (IT) 0);
   
    for (IT k = 0 ; k < nnz ; ++k)
    {
        IT tmp =  triples[k].col;
        work [ tmp ]++ ;		// column counts (i.e, w holds the "col difference array")
	}

	if(nnz > 0)
	{
		colptr[cols] = CumulativeSum (work, cols) ;		// cumulative sum of w
        	copy(work, work+cols, colptr);
		IT last;
		for (IT k = 0 ; k < nnz ; ++k)
        	{
			tosort[ work[triples[k].col]++] = make_pair( triples[k].row, triples[k].val);
		}
		#pragma omp parallel for
		for(int i=0; i< cols; ++i)
		{
			sort(tosort.begin() + colptr[i], tosort.begin() + colptr[i+1]);

			typename vector<pair<IT,NT> >::iterator itr;	// iterator is a dependent name
			IT ind;
			for(itr = tosort.begin() + colptr[i], ind = colptr[i]; itr != tosort.begin() + colptr[i+1]; ++itr, ++ind)
			{
				rowids[ind] = itr->first;
				values[ind] = itr->second;
			}
        	}
	}
	delete [] work;
}


template <class IT, class NT>
template <typename AddOperation>
void CSC<IT,NT>::MergeDuplicates (AddOperation addop)
{
    vector<IT> diff(cols,0);
    std::adjacent_difference(colptr+1, colptr+cols+1, diff.begin());
    
    vector<vector<IT>> v_rowids;
    vector<vector<NT>> v_values;

    if(nnz > 0)
    {
        //#pragma omp parallel for
        for(int i=0; i<cols; ++i)
        {
            for(size_t j=colptr[i]; j<colptr[i+1]; ++j)
            {
                v_rowids[i].push_back(rowids[j]); // SEG FAULT
                v_values[i].push_back(values[j]);
                while(j < colptr[i+1]-1 && rowids[j] == rowids[j+1])
                {
                    v_values[i].back() = addop(v_values[i].back(), values[j+1]);
                    j++;    // increment j
                    diff[i]--;
                }
            }
        }
    }

    colptr[cols] = CumulativeSum (diff.data(), cols) ;	// cumulative sum of diff
    copy(diff.begin(), diff.end(), colptr);  // update the column pointers
    delete [] rowids;
    delete [] values;
    cout << "old number of nonzeros before merging : " << nnz << endl; 
    nnz  = colptr[cols];
    cout << "new number of nonzeros after merging : " << nnz << endl; 
    
    rowids = new IT[nnz];
    values = new NT[nnz];
    
//#pragma omp parallel for
    for(int i=0; i<cols; ++i)
    {
        copy(v_rowids[i].begin(), v_rowids[i].end(), rowids+colptr[i]);
        copy(v_values[i].begin(), v_values[i].end(), values+colptr[i]);
    }
}
//! this version handles duplicates in the input
template <class IT, class NT>
template <typename AddOperation>
CSC<IT,NT>::CSC (vector< tuple<IT,IT,NT> > & tuple, IT m, IT n, AddOperation addop): rows(m), cols(n)
{
    nnz = tuple.size(); // there might be duplicates
    // cout << nnz << endl; 281528
    // cout << cols << endl; 1636281
    colptr = new IT[cols]();
    rowids = new IT[nnz];
    values = new NT[nnz];
    vector< pair<IT,NT> > tosort (nnz);

    IT *work = new IT[cols](); // workspace
    std::fill(work, work+cols, (IT) 0); 

    for (IT k = 0; k < nnz; ++k)
    {
        IT tmp = get<1>(tuple[k]); 
        work[tmp]++;	        	// column counts (i.e, work holds the "col difference array") 
    }
    if(nnz > 0)
    {
        colptr[cols] = CumulativeSum (work, cols);		// cumulative sum of work
        copy(work, work+cols, colptr);
        IT last;
        for (IT k = 0 ; k < nnz ; ++k)
        {
            tosort[work[get<1>(tuple[k])]++] = make_pair(get<0>(tuple[k]), get<2>(tuple[k]));
        }
//#pragma omp parallel for
        for(int i=0; i< cols; ++i)
        {
            sort(tosort.begin() + colptr[i], tosort.begin() + colptr[i+1]);
            typename vector<pair<IT,NT> >::iterator itr;	// iterator is a dependent name
            IT ind;
            for(itr = tosort.begin() + colptr[i], ind = colptr[i]; itr != tosort.begin() + colptr[i+1]; ++itr, ++ind)
            {
                rowids[ind] = itr->first;
                values[ind] = itr->second;
            }
        }
    }

    delete [] work;
    cerr << "Before MergeDuplicates()" << endl;
    MergeDuplicates(addop);
}


// Construct a Csc object from parallel arrays
template <class IT, class NT>
CSC<IT,NT>::CSC(IT * ri, IT * ci, NT * val, IT mynnz, IT m, IT n):nnz(mynnz),rows(m),cols(n)
{
    assert(nnz != 0 && rows != 0);
    colptr = new IT[cols+1];
    rowids = new IT[nnz];
    values = new NT[nnz];
    vector< pair<IT,NT> > tosort (nnz); 

    IT * work = new IT[cols];	// workspace
    std::fill(work, work+cols, (IT) 0);

    for (IT k = 0; k < nnz; ++k)
    {
        IT tmp = ci[k];
        work[ tmp ]++;		// column counts (i.e, w holds the "col difference array")
    }
    if(nnz > 0)
    {
	colptr[cols] = CumulativeSum (work, cols) ;		// cumulative sum of w
        copy(work, work+cols, colptr);
	IT last;
	for (IT k = 0 ; k < nnz ; ++k)
        {
		tosort[ work[ci[k]]++] = make_pair( ri[k], val[k]);
	}
	#pragma omp parallel for
	for(int i=0; i< cols; ++i)
	{
		sort(tosort.begin() + colptr[i], tosort.begin() + colptr[i+1]);
		typename vector<pair<IT,NT> >::iterator itr;	// iterator is a dependent name
		IT ind;
		for(itr = tosort.begin() + colptr[i], ind = colptr[i]; itr != tosort.begin() + colptr[i+1]; ++itr, ++ind)
		{
			rowids[ind] = itr->first;
			values[ind] = itr->second;
		}
        }
    }
    delete [] work;
}

// check if sorted within columns
template <class IT, class NT>
void CSC<IT,NT>::Sorted()
{
	// CSR<IT,NT> csr (*this);    	// convert to csr
	// this = new CSC<IT,NT>(csr);	// convert back to csc
	bool sorted = true;
	for(IT i=0; i< cols; ++i)
	{
		sorted &= my_is_sorted (rowids + colptr[i], rowids + colptr[i+1], std::less<IT>());
	}
	cout << "Sorted ? " << sorted << endl;
}

template <class IT, class NT>
bool CSC<IT,NT>::operator==(const CSC<IT,NT> & rhs)
{
  if(nnz != rhs.nnz || rows  != rhs.rows || cols != rhs.cols) return false;
  bool same = std::equal(colptr, colptr+cols+1, rhs.colptr); 
  same = same && std::equal(rowids, rowids+nnz, rhs.rowids);
  
   bool samebefore = same;
 ErrorTolerantEqual<NT> epsilonequal(EPSILON);
  same = same && std::equal(values, values+nnz, rhs.values, epsilonequal );
  if(samebefore && (!same))
  {
#ifdef DEBUG
	vector<NT> error(nnz);
	transform(values, values+nnz, rhs.values, error.begin(), absdiff<NT>());
	vector< pair<NT, NT> > error_original_pair(nnz);
	for(IT i=0; i < nnz; ++i)
		error_original_pair[i] = make_pair(error[i], values[i]);
	if(error_original_pair.size() > 10)	// otherwise would crush for small data
	{
		partial_sort(error_original_pair.begin(), error_original_pair.begin()+10, error_original_pair.end(), greater< pair<NT,NT> >());
		cout << "Highest 10 different entries are: " << endl;
		for(IT i=0; i < 10; ++i)
			cout << "Diff: " << error_original_pair[i].first << " on " << error_original_pair[i].second << endl;
	}
	else
	{
		sort(error_original_pair.begin(), error_original_pair.end(), greater< pair<NT,NT> >());
		cout << "Highest different entries are: " << endl;
		for(typename vector< pair<NT, NT> >::iterator it=error_original_pair.begin(); it != error_original_pair.end(); ++it)
			cout << "Diff: " << it->first << " on " << it->second << endl;
	}
#endif
}
  return same;
}

template <class IT, class NT>
void CSC<IT,NT>::intersect (const IT* rowids_in, const NT* values_in, const IT len_in,
                            const IT* ri, const IT len_ri,
                            IT* rowids_out, NT* values_out, IT* len_out)
{
    IT maxlen = len_in>len_ri ? len_in : len_ri;
    double r = len_in>len_ri ? (double)len_in/len_ri : (double)len_ri/len_in;
    //if(log2(maxlen) < r) // linear scan is asymptotically better
    {
        IT idx=0;
        for(int j=0, k=0; j<len_in && k < len_ri;)
        {
            if(ri[k] < rowids_in[j]) k++;
            else if(ri[k] > rowids_in[j]) j++;
            else //(ri[k]==rowids[j])
            {
                values_out[idx] = values_in[j];
                rowids_out[idx++] = rowids_in[j];
                k++;
                j++; // repeated rows are not allowed
            }
        }
        *len_out = idx;
    
    }
    //else // use finger search
    {
    }
    
}

template <class IT, class NT>
CSC<IT,NT> CSC<IT,NT>::SpRef2 (const IT* ri, const IT rilen, const IT* ci, const IT cilen)
{
    if( cilen>0 && ci[cilen-1] > cols)
    {
        cerr << "Col indices out of bounds" << endl;
        abort();
    }
    if( rilen>0 && ri[rilen-1] > rows)
    {
        cerr << "Row indices out of bounds" << endl;
        abort();
    }
    
    
    // count nnz(A[,:J])
    IT nnz_ci = 0;
    for(int i=0; i<cilen; i++)
    {
        nnz_ci = nnz_ci + colptr[ci[i]+1] - colptr[ci[i]];

    }
    


    IT* rowids_out = new IT[nnz_ci];
    NT* values_out = new NT[nnz_ci];
    IT* len_out = new IT[cilen];
    IT idx=0;
    for(int i=0; i<cilen; i++)
    {
        IT cidx1 = colptr[ci[i]];
        IT cidx2 = colptr[ci[i]+1];
        
        intersect (&rowids[cidx1], &values[cidx1], cidx2 - cidx1,ri, rilen,
                   &rowids_out[cidx1], &values_out[cidx1], &len_out[i]);
    }
    
    
    CSC C;
    C.rows = rilen;
    C.cols = cilen;
    C.colptr = new IT[C.cols+1];
    C.colptr[0] = 0;
        
    for(int i=0; i < C.cols; ++i)
    {
        C.colptr[i+1] = C.colptr[i] + len_out[i];
    }
    C.nnz = C.colptr[C.cols];
    C.rowids = new IT[C.nnz];
    C.values = new NT[C.nnz];
        

    for(int i=0; i< C.cols; ++i)         // combine step
    {
        IT cidx1 = colptr[ci[i]];
        IT cidx2 = cidx1 + len_out[i];
        copy(&rowids_out[cidx1], &rowids_out[cidx2], C.rowids + C.colptr[i]);
        copy(&values_out[cidx1], &values_out[cidx2], C.values + C.colptr[i]);
    }
    return C;
}


// write genereal purpose set-intersect
// binary search is faster is one of the vectors is very large


// we assume that ri and ci are sorted in ascending order
// also assume that matrix sorted within column
// output is another CSC
// note that ri and ci might have repeated entries
// behaviour is exactly similar to the matlab implementation
template <class IT, class NT>
CSC<IT,NT> CSC<IT,NT>::SpRef (const vector<IT> & ri, const vector<IT> & ci)
{
    if( (!ci.empty()) && (ci.back() > cols))
    {
        cerr << "Col indices out of bounds" << endl;
        abort();
    }
    if( (!ri.empty()) && (ri.back() > rows))
    {
        cerr << "Row indices out of bounds" << endl;
        abort();
    }
    
   
    // first, count nnz in the result matrix
    IT refnnz = 0;
    for(int i=0; i<ci.size(); i++)
    {
        IT j = colptr[ci[i]], k=0;
        IT endIdx = colptr[ci[i]+1];
        while(j<endIdx && k < ri.size())
        {
            //cout << j << "=" << rowids[j] << " :: " << k << "=" << ri[k] << " \n";
            if(ri[k]<rowids[j]) k++;
            else if(ri[k]>rowids[j]) j++;
            else //(ri[k]==rowids[j])
            {
                
                refnnz++;
                k++;
                //j++;  // wait for the next iteration of the inner loop to alow reapted rows
            }
        }
    }


    // Next, allocate memory and save the result matrix
    // This two-step implementation is better for multithreading
    CSC refmat(refnnz, ri.size(), ci.size(), 0);
    refmat.colptr[0] = 0;
    IT idx=0;
    for(int i=0; i<ci.size(); i++)
    {
        IT j = colptr[ci[i]], k=0;
        IT endIdx = colptr[ci[i]+1];
        while(j<endIdx && k < ri.size())
        {
            if(ri[k]<rowids[j]) k++;
            else if(ri[k]>rowids[j]) j++;
            else //(ri[k]==rowids[j])
            {
                refmat.values[idx] = values[j];
                refmat.rowids[idx++] = rowids[j];
                k++;
                //j++; // wait for the next iteration of the inner loop to alow reapted rows
            }
        }
        refmat.colptr[i+1] = idx;
    }
    
    return refmat;
}


// write genereal purpose set-intersect
// binary search is faster is one of the vectors is very large


// we assume that ri and ci are sorted in ascending order
// also assume that matrix sorted within column
// output is another CSC
// note that ri and ci might have repeated entries
// behaviour is exactly similar to the matlab implementation
template <class IT, class NT>
CSC<IT,NT> CSC<IT,NT>::SpRef1 (const vector<IT> & ri, const vector<IT> & ci)
{
    if( (!ci.empty()) && (ci.back() > cols))
    {
        cerr << "Col indices out of bounds" << endl;
        abort();
    }
    if( (!ri.empty()) && (ri.back() > rows))
    {
        cerr << "Row indices out of bounds" << endl;
        abort();
    }
    
    BitMap bmap(ri.size()); // space requirement n bits
    bmap.reset(); // this is time consuming .....
    for(int i=0; i<ri.size(); i++)
    {
        bmap.set_bit(ri[i]);
    }
    
    // first, count nnz in the result matrix
    IT refnnz = 0;
    for(int i=0; i<ci.size(); i++)
    {
        IT endIdx = colptr[ci[i]+1];
        for(IT j=colptr[ci[i]]; j<endIdx; j++)
        {
            if(bmap.get_bit(rowids[j])) refnnz++;
        }
    }
    
    
    // Next, allocate memory and save the result matrix
    // This two-step implementation is better for multithreading
    CSC refmat(refnnz, ri.size(), ci.size(), 0);
    refmat.colptr[0] = 0;
    IT idx=0;
    for(int i=0; i<ci.size(); i++)
    {
        IT endIdx = colptr[ci[i]+1];
        for(IT j=colptr[ci[i]]; j<endIdx; j++)
        {
            if(bmap.get_bit(rowids[j]))
            {
                refmat.values[idx] = values[j];
                refmat.rowids[idx++] = rowids[j];
            }
        }
        refmat.colptr[i+1] = idx;
    }
    
    
    return refmat;
}