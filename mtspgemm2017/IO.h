#ifndef _IO_SPGEMM_H
#define _IO_SPGEMM_H

#include "Triple.h"
#include "CSC.h"

#define READBUFFER (512 * 1024 * 1024)  // in MB

template <typename IT, typename NT>
int ReadBinary(string filename, CSC<IT,NT> * & csc)
{
    FILE * f = fopen(filename.c_str(), "r");
    if(!f)
    {
        cerr << "Problem reading binary input file" << filename << endl;
        return -1;
    }
    IT m,n,nnz;
    fread(&m, sizeof(IT), 1, f);
    fread(&n, sizeof(IT), 1, f);
    fread(&nnz, sizeof(IT), 1, f);
    
    if (m <= 0 || n <= 0 || nnz <= 0)
    {
        cerr << "Problem with matrix size in binary input file" << filename << endl;
        return -1;
    }
    double start = omp_get_wtime( );
    cout << "Reading matrix with dimensions: "<< m << "-by-" << n <<" having "<< nnz << " nonzeros" << endl;
    
    IT * rowindices = new IT[nnz];
    IT * colindices = new IT[nnz];
    NT * vals = new NT[nnz];
    
    size_t rows = fread(rowindices, sizeof(IT), nnz, f);
    size_t cols = fread(colindices, sizeof(IT), nnz, f);
    size_t nums = fread(vals, sizeof(NT), nnz, f);
    
    if(rows != nnz || cols != nnz || nums != nnz)
    {
        cerr << "Problem with FREAD, aborting... " << endl;
        return -1;
    }
    double end = omp_get_wtime( );
    printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

    fclose(f);
    
    csc = new CSC<IT,NT>(rowindices, colindices, vals , nnz, m, n);
    
    delete [] rowindices;
    delete [] colindices;
    delete [] vals;
    return 1;
}

template <typename IT, typename NT>
int ReadASCII(string filename, CSC<IT,NT> * & csc)
{
    double start = omp_get_wtime( );
    ifstream infile(filename.c_str());
    char line[256];
    char c = infile.get();
    while(c == '%')
    {
        infile.getline(line,256);
        c = infile.get();
    }
    infile.unget();
    IT m,n,nnz;
    infile >> m >> n >> nnz;	// #{rows}-#{cols}-#{nonzeros}
    //cout << m << " " << n << " " << nnz << endl;
    
    Triple<IT,NT> * triples = new Triple<IT,NT>[nnz];
    if (infile.is_open())
    {
        IT cnz = 0;	// current number of nonzeros
        while (! infile.eof() && cnz < nnz)
        {
            infile >> triples[cnz].row >> triples[cnz].col >> triples[cnz].val;	// row-col-value
            triples[cnz].row--;
            triples[cnz].col--;
            ++cnz;
        }
        assert(cnz == nnz);
    }
    
    double end = omp_get_wtime( );
    printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
	
    cout << "converting to csc ... " << endl;
    csc= new CSC<IT,NT>(triples, nnz, m, n);
    csc->totalcols = n;
    delete [] triples;
    return 1;
}



#endif
