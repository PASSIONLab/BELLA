
template <typename KIND, typename CIND>  // kmer_index, count_index
void WriteToDisk(const vector<vector<tuple<KIND, KIND, CIND>>> & alltranstuples,
                 const CuckooDict<KIND> & countsreliable,
                 KIND readcount, KIND tuplecount)
{
	cout << "Writing to disk" << endl;
	string filename = "readbykmers.mtx";
	vector<std::stringstream> ss(MAXTHREADS);
	vector<std::string> text(MAXTHREADS);
	vector<int64_t> bytes(MAXTHREADS);

	ss[0] << readcount << '\t' << countsreliable.size() << '\t' << tuplecount << '\n';  
	#pragma omp parallel for
	for(int t=0; t<MAXTHREADS; ++t)
        {	
		for( auto tlp: alltranstuples[t])
		{
			ss[t] << get<1>(tlp)+1 << '\t' <<  get<0>(tlp)+1 <<  '\t' << get<2>(tlp) << '\n';
		}
		text[t] = ss[t].str();
		bytes[t] = text[t].size();
		ss[t].clear();
	}
	std::ofstream ofs(filename.c_str(), std::ios::binary | std::ios::out);
	vector<int64_t> bytesuntil(MAXTHREADS+1, 0);
	std::partial_sum(bytes.begin(), bytes.end(), bytesuntil.begin()+1);	
	ofs.seekp(bytes[MAXTHREADS-1] - 1);
	ofs.write("", 1);   // this will likely create a sparse file so the actual disks won't spin yet
	ofs.close();
	struct stat st;     // get file size
        if (stat(filename.c_str(), &st) != -1)
        {
		std::cout << "File is actually " << st.st_size << " bytes" << endl;
        }

	#pragma omp parallel for
        for(int t=0; t<MAXTHREADS; ++t)
        {
		FILE *ffinal = fopen(filename.c_str(), "rb+");
		fseek (ffinal , bytesuntil[t] , SEEK_SET );
     		fwrite(text[t].c_str(),1, bytes[t] ,ffinal);
     		fflush(ffinal);
     		fclose(ffinal);
	}
	cout << "Output of k-mer x read matrix written" << endl;
}
