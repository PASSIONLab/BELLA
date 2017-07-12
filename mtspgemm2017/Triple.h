#ifndef _TRIPLE_H_
#define _TRIPLE_H_

#include <iostream>
#include <functional>
using namespace std;

template <class IT, class NT>
struct Triple
{
	IT row;	//  row index
	IT col;	//  col index
	NT val;		//  value
	Triple(IT myrow, IT mycol, NT myval):row(myrow),col(mycol),val(myval) {};
	Triple():row(0),col(0),val(0) {};


	class ScalarReadSaveHandler
	{
	public:
		template <typename c, typename t>
		void save(std::basic_ostream<c,t>& os, const NT& v)
		{
			os << v;
		}
	};

	void ParallelWrite(const string & filename, bool onebased, bool includeindices = true) { ParallelWrite(filename, onebased, ScalarReadSaveHandler(), includeindices); };	
	
	template <class HANDLER>	
	void ParallelWrite(const string & filename, bool onebased, HANDLER handler, bool includeindices)
	{
		#pragma omp parallel
   		{
        		long int fpos, end_fpos;
        		int myrank = omp_get_thread_num();
        		int nprocs = omp_get_num_threads();

       			IT totalLength = TotalLength();
			IT totalNNZ = getnnz();

			std::stringstream ss;
			if(myrank == 0)
			{
				ss << totalLength << '\t' << totalNNZ << '\n';	// rank-0 has the header
			}	
			IT entries =  getlocnnz();
			IT sizeuntil = 0;
			MPI_Exscan( &entries, &sizeuntil, 1, MPIType<IT>(), MPI_SUM, commGrid->GetWorld() );
	if(myrank == 0) sizeuntil = 0;	// because MPI_Exscan says the recvbuf in process 0 is undefined

	if(includeindices)
	{
		if(onebased)	sizeuntil += 1;	// increment by 1
		
		for(IT i=0; i< entries; ++i)
		{
			ss << ind[i]+sizeuntil << '\t';
			handler.save(ss, num[i]);
			ss << '\n';
		}
	}
	else	// the base doesn't matter if we don't include indices
	{
		IT dummy = 0;	// dummy because we don't want indices to be printed
		for(IT i=0; i< entries; ++i)
		{
			handler.save(ss, num[i], dummy);
			ss << '\n';
		}
	}
	
	std::string text = ss.str();

	int64_t * bytes = new int64_t[nprocs];
    	bytes[myrank] = text.size();
    	MPI_Allgather(MPI_IN_PLACE, 1, MPIType<int64_t>(), bytes, 1, MPIType<int64_t>(), commGrid->GetWorld());
	int64_t bytesuntil = accumulate(bytes, bytes+myrank, static_cast<int64_t>(0));
	int64_t bytestotal = accumulate(bytes, bytes+nprocs, static_cast<int64_t>(0));

    	if(myrank == 0)	// only leader rights the original file with no content
    	{
		std::ofstream ofs(filename.c_str(), std::ios::binary | std::ios::out);
		cout << "Creating file with " << bytestotal << " bytes" << endl;
    		ofs.seekp(bytestotal - 1);
    		ofs.write("", 1);	// this will likely create a sparse file so the actual disks won't spin yet
		ofs.close();
   	}
     	MPI_Barrier(commGrid->GetWorld());
    
	struct stat st;     // get file size
    	if (stat(filename.c_str(), &st) == -1)
    	{
       		MPI_Abort(commGrid->GetWorld(), NOFILE);
    	}
	if(myrank == nprocs-1)	// let some other processor do the testing
	{
		cout << "File is actually " << st.st_size << " bytes seen from process " << myrank << endl;	
	}

    	FILE *ffinal;
	if ((ffinal = fopen(filename.c_str(), "rb+")) == NULL)	// then everyone fills it
        {
		printf("COMBBLAS: Vector output file %s failed to open at process %d\n", filename.c_str(), myrank);
            	MPI_Abort(commGrid->GetWorld(), NOFILE);
       	}
	fseek (ffinal , bytesuntil , SEEK_SET );
	fwrite(text.c_str(),1, bytes[myrank] ,ffinal);
	fflush(ffinal);
	fclose(ffinal);
	delete [] bytes;
};

#endif
