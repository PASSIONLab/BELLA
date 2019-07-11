//==================================================================
// Title:  Cuda x-drop seed-and-extend alignment algorithm
// Author: G. Guidi, A. Zeni
// Date:   6 March 2019
//==================================================================
#define N_THREADS 512 
#define N_BLOCKS 500000 
#define MIN -32768
#define BYTES_INT 4
#define XDROP 21
// #define N_STREAMS 60
#define MAX_SIZE_ANTIDIAG 8000
#define MAX_GPUS 8

//trying to see if the scoring scheme is a bottleneck in some way
#define MATCH     1
#define MISMATCH -1
#define GAP_EXT  -1
#define GAP_OPEN -1
#define UNDEF -32767
#define NOW std::chrono::high_resolution_clock::now()

using namespace std;
using namespace chrono;

#define cudaErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char* file, int line, bool abort = true){

	if(code != cudaSuccess){
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
	if(abort) exit(code);
	}
}

enum ExtensionDirectionL
{
	EXTEND_NONEL  = 0,
	EXTEND_LEFTL  = 1,
	EXTEND_RIGHTL = 2,
	EXTEND_BOTHL  = 3
};

__inline__ __device__ void warpReduce(volatile short *input,
										  int myTId){
		input[myTId] = (input[myTId] > input[myTId + 32]) ? input[myTId] : input[myTId + 32]; 
		input[myTId] = (input[myTId] > input[myTId + 16]) ? input[myTId] : input[myTId + 16];
		input[myTId] = (input[myTId] > input[myTId + 8]) ? input[myTId] : input[myTId + 8]; 
		input[myTId] = (input[myTId] > input[myTId + 4]) ? input[myTId] : input[myTId + 4];
		input[myTId] = (input[myTId] > input[myTId + 2]) ? input[myTId] : input[myTId + 2];
		input[myTId] = (input[myTId] > input[myTId + 1]) ? input[myTId] : input[myTId + 1];
}

__inline__ __device__ short reduce_max(short *input, int dim){
	unsigned int myTId = threadIdx.x;   
	if(dim>32){
		for(int i = N_THREADS/2; i >32; i>>=1){
			if(myTId < i){
						input[myTId] = (input[myTId] > input[myTId + i]) ? input[myTId] : input[myTId + i];
			}__syncthreads();
		}//__syncthreads();
	}
	if(myTId<32)
		warpReduce(input, myTId);
	__syncthreads();
	return input[0];
}

__inline__ __device__ void updateExtendedSeedL(SeedL &seed,
					ExtensionDirectionL direction, //as there are only 4 directions we may consider even smaller data types
					int cols,
					int rows,
					int lowerDiag,
					int upperDiag)
{
	if (direction == EXTEND_LEFTL)
	{
		int beginDiag = seed.beginDiagonal;
		// Set lower and upper diagonals.
		
		if (getLowerDiagonal(seed) > beginDiag + lowerDiag)
			setLowerDiagonal(seed, beginDiag + lowerDiag);
		if (getUpperDiagonal(seed) < beginDiag + upperDiag)
			setUpperDiagonal(seed, beginDiag + upperDiag);

		// Set new start position of seed.
		seed.beginPositionH -= rows;
		seed.beginPositionV -= cols;
	} else {  // direction == EXTEND_RIGHTL
		// Set new lower and upper diagonals.
		int endDiag = seed.endDiagonal;
		if (getUpperDiagonal(seed) < endDiag - lowerDiag)
			setUpperDiagonal(seed, (endDiag - lowerDiag));
		if (getLowerDiagonal(seed) > (endDiag - upperDiag))
			setLowerDiagonal(seed, endDiag - upperDiag);

		// Set new end position of seed.
		seed.endPositionH += rows;
		seed.endPositionV += cols;
		
	}
}

__inline__ __device__ void computeAntidiag(short *antiDiag1,
									short *antiDiag2,
									short *antiDiag3,
									char* querySeg,
									char* databaseSeg,
									int &best,
									int &scoreDropOff,
									int &cols,
									int &rows,
									int &minCol,
									int &maxCol,
									int &antiDiagNo,
									int &offset1,
									int &offset2,
									ExtensionDirectionL direction
									){
	int tid = threadIdx.x;
	
	for(int i = 0; i < maxCol; i+=N_THREADS){

		int col = tid + minCol + i;
		int queryPos, dbPos;
		
		queryPos = col - 1;
		dbPos = col + rows - antiDiagNo - 1;
		
		/*if(direction == EXTEND_LEFTL){
			queryPos = cols - 1 - col;
			dbPos = rows - 1 + col - antiDiagNo;
		}else{//EXTEND RIGHT
			queryPos = col - 1;
			dbPos = antiDiagNo - col - 1;
		}*/

		if(col < maxCol){
		
			int tmp = max_logan(antiDiag2[col-offset2],antiDiag2[col-offset2-1]) + GAP_EXT;
		
			int score = (querySeg[queryPos] == databaseSeg[dbPos]) ? MATCH : MISMATCH;
			
			tmp = max_logan(antiDiag1[col-offset1-1]+score,tmp);
			
			antiDiag3[tid+1+i] = (tmp < best - scoreDropOff) ? UNDEF : tmp;
		
		}
	}
}

__inline__ __device__ void calcExtendedLowerDiag(int &lowerDiag,
			  int const &minCol,
			  int const &antiDiagNo)
{
	int minRow = antiDiagNo - minCol;
	if (minCol - minRow < lowerDiag)
		lowerDiag = minCol - minRow;
}

__inline__ __device__ void calcExtendedUpperDiag(int &upperDiag,
			  int const &maxCol,
			  int const &antiDiagNo)
{
	int maxRow = antiDiagNo + 1 - maxCol;
	if (maxCol - 1 - maxRow > upperDiag)
		upperDiag = maxCol - 1 - maxRow;
}

__inline__ __device__ void initAntiDiag3(short *antiDiag3,
						   int &a3size,
						   int const &offset,
						   int const &maxCol,
						   int const &antiDiagNo,
						   int const &minScore,
						   int const &gapCost,
						   int const &undefined)
{
	a3size = maxCol + 1 - offset;

	antiDiag3[0] = undefined;
	antiDiag3[maxCol - offset] = undefined;

	if (antiDiagNo * gapCost > minScore)
	{
		if (offset == 0) // init first column
			antiDiag3[0] = antiDiagNo * gapCost;
		if (antiDiagNo - maxCol == 0) // init first row
			antiDiag3[maxCol - offset] = antiDiagNo * gapCost;
	}
	//return offset;
}

__inline__ __device__ void initAntiDiags(
			   short *antiDiag1,
			   short *antiDiag2,
			   short *antiDiag3,
			   int &a2size,
			   int &a3size,
			   int const &dropOff,
			   int const &gapCost,
			   int const &undefined)
{
	//antiDiag2.resize(1);
	a2size = 1;

	//resize(antiDiag2, 1);
	antiDiag2[0] = 0;

	//antiDiag3.resize(2);
	a3size = 2;

	//if (-gapCost > dropOff)
	//{
	// 	antiDiag3[0] = undefined;
	// 	antiDiag3[1] = undefined;
	//}
	//else
	//{
		antiDiag3[0] = gapCost;
		antiDiag3[1] = gapCost;
	//}
}

__global__ void extendSeedLGappedXDropOneDirectionShared(
		SeedL *seed,
		char *querySegArray,
		char *databaseSegArray,
		ExtensionDirectionL direction,
		int scoreDropOff,
		int *res,
		int *offsetQuery,
		int *offsetTarget,
		int offAntidiag
		)
{
	int myId = blockIdx.x;
	int myTId = threadIdx.x;
	char *querySeg;
	char *databaseSeg;

	if(myId==0){
		querySeg = querySegArray;
		databaseSeg = databaseSegArray;
	}
	else{
		querySeg = querySegArray + offsetQuery[myId-1];
		databaseSeg = databaseSegArray + offsetTarget[myId-1];
	}

	extern __shared__ short antidiagonals[]; //decomment this for shared/ comment for global
	short* antiDiag1 = &antidiagonals[0]; //decomment this for shared/ comment for global
	short* antiDiag2 = &antiDiag1[offAntidiag];
	short* antiDiag3 = &antiDiag2[offAntidiag];


	SeedL mySeed(seed[myId]);	
	//dimension of the antidiagonals
	int a1size = 0, a2size = 0, a3size = 0;
	int cols, rows;

	if(myId == 0){
			cols = offsetQuery[myId]+1;
			rows = offsetTarget[myId]+1;
	}
	else{
			cols = offsetQuery[myId]-offsetQuery[myId-1]+1;
			rows = offsetTarget[myId]-offsetTarget[myId-1]+1;
	}

	if (rows == 1 || cols == 1)
		return;

	//printf("%d\n", gapCost);
	//int undefined = UNDEF;
	int minCol = 1;
	int maxCol = 2;

	int offset1 = 0; // number of leading columns that need not be calculated in antiDiag1
	int offset2 = 0; //                                                       in antiDiag2
	int offset3 = 0; //                                                       in antiDiag3

	initAntiDiags(antiDiag1,antiDiag2, antiDiag3, a2size, a3size, scoreDropOff, GAP_EXT, UNDEF);
	int antiDiagNo = 1; // the currently calculated anti-diagonal

	int best = 0; // maximal score value in the DP matrix (for drop-off calculation)

	int lowerDiag = 0;
	int upperDiag = 0;

	while (minCol < maxCol)
	{	

		
		++antiDiagNo;

		//antidiagswap
		//antiDiag2 -> antiDiag1
		//antiDiag3 -> antiDiag2
		//antiDiag1 -> antiDiag3
		short *t = antiDiag1;
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag3 = t;
		int t_l = a1size;
		a1size = a2size;
		a2size = a3size;
		a3size = t_l;
		offset1 = offset2;
		offset2 = offset3;
		offset3 = minCol-1;
		__shared__ short temp[N_THREADS];
		initAntiDiag3(antiDiag3, a3size, offset3, maxCol, antiDiagNo, best - scoreDropOff, GAP_EXT, UNDEF);
		
		computeAntidiag(antiDiag1, antiDiag2, antiDiag3, querySeg, databaseSeg, best, scoreDropOff, cols, rows, minCol, maxCol, antiDiagNo, offset1, offset2, direction);	 	
		__syncthreads();	
	
		int tmp, antiDiagBest = UNDEF;	
		for(int i=0; i<a3size; i+=N_THREADS){
			int size = a3size-i;
			
			if(myTId<N_THREADS){
				temp[myTId] = (myTId<size) ? antiDiag3[myTId+i]:UNDEF;				
			}
			__syncthreads();
			
			tmp = reduce_max(temp,size);
			antiDiagBest = (tmp>antiDiagBest) ? tmp:antiDiagBest;

		}
		best = (best > antiDiagBest) ? best : antiDiagBest;
		
		while (minCol - offset3 < a3size && antiDiag3[minCol - offset3] == UNDEF &&
			   minCol - offset2 - 1 < a2size && antiDiag2[minCol - offset2 - 1] == UNDEF)
		{
			++minCol;
		}

		// Calculate new maxCol
		while (maxCol - offset3 > 0 && (antiDiag3[maxCol - offset3 - 1] == UNDEF) &&
									   (antiDiag2[maxCol - offset2 - 1] == UNDEF))
		{
			--maxCol;
		}
		++maxCol;

		// Calculate new lowerDiag and upperDiag of extended seed
		calcExtendedLowerDiag(lowerDiag, minCol, antiDiagNo);
		calcExtendedUpperDiag(upperDiag, maxCol - 1, antiDiagNo);
		
		// end of databaseSeg reached?
		minCol = (minCol > (antiDiagNo + 2 - rows)) ? minCol : (antiDiagNo + 2 - rows);
		// end of querySeg reached?
		maxCol = (maxCol < cols) ? maxCol : cols;
		//}
	}

	int longestExtensionCol = a3size + offset3 - 2;
	int longestExtensionRow = antiDiagNo - longestExtensionCol;
	int longestExtensionScore = antiDiag3[longestExtensionCol - offset3];
	
	if (longestExtensionScore == UNDEF)
	{
		if (antiDiag2[a2size -2] != UNDEF)
		{
			// reached end of query segment
			longestExtensionCol = a2size + offset2 - 2;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
			
		}
		else if (a2size > 2 && antiDiag2[a2size-3] != UNDEF)
		{
			// reached end of database segment
			longestExtensionCol = a2size + offset2 - 3;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
			
		}
	}


	if (longestExtensionScore == UNDEF){

		// general case
		for (int i = 0; i < a1size; ++i){

			if (antiDiag1[i] > longestExtensionScore){

				longestExtensionScore = antiDiag1[i];
				longestExtensionCol = i + offset1;
				longestExtensionRow = antiDiagNo - 2 - longestExtensionCol;
			
			}
		}
	}
	
	if (longestExtensionScore != UNDEF)
		updateExtendedSeedL(mySeed, direction, longestExtensionCol, longestExtensionRow, lowerDiag, upperDiag);
	seed[myId] = mySeed;
	res[myId] = longestExtensionScore;

}

__global__ void extendSeedLGappedXDropOneDirectionGlobal(
		SeedL *seed,
		char *querySegArray,
		char *databaseSegArray,
		ExtensionDirectionL direction,
		int scoreDropOff,
		int *res,
		int *offsetQuery,
		int *offsetTarget,
		int offAntidiag,
		short *antidiag
		)
{
	int myId = blockIdx.x;
	int myTId = threadIdx.x;
	char *querySeg;
	char *databaseSeg;

	if(myId==0){
		querySeg = querySegArray;
		databaseSeg = databaseSegArray;
	}
	else{
		querySeg = querySegArray + offsetQuery[myId-1];
		databaseSeg = databaseSegArray + offsetTarget[myId-1];
	}

	short *antiDiag1 = &antidiag[myId*offAntidiag*3]; 
	short* antiDiag2 = &antiDiag1[offAntidiag];
		short* antiDiag3 = &antiDiag2[offAntidiag];


	SeedL mySeed(seed[myId]);	
	//dimension of the antidiagonals
	int a1size = 0, a2size = 0, a3size = 0;
	int cols, rows;

	if(myId == 0){
			cols = offsetQuery[myId]+1;
			rows = offsetTarget[myId]+1;
	}
	else{
			cols = offsetQuery[myId]-offsetQuery[myId-1]+1;
			rows = offsetTarget[myId]-offsetTarget[myId-1]+1;
	}

	if (rows == 1 || cols == 1)
		return;

	//printf("%d\n", gapCost);
	//int undefined = UNDEF;
	int minCol = 1;
	int maxCol = 2;

	int offset1 = 0; // number of leading columns that need not be calculated in antiDiag1
	int offset2 = 0; //                                                       in antiDiag2
	int offset3 = 0; //                                                       in antiDiag3

	initAntiDiags(antiDiag1,antiDiag2, antiDiag3, a2size, a3size, scoreDropOff, GAP_EXT, UNDEF);
	int antiDiagNo = 1; // the currently calculated anti-diagonal

	int best = 0; // maximal score value in the DP matrix (for drop-off calculation)

	int lowerDiag = 0;
	int upperDiag = 0;

	while (minCol < maxCol)
	{	

		
		++antiDiagNo;

		//antidiagswap
		//antiDiag2 -> antiDiag1
		//antiDiag3 -> antiDiag2
		//antiDiag1 -> antiDiag3
		short *t = antiDiag1;
		antiDiag1 = antiDiag2;
		antiDiag2 = antiDiag3;
		antiDiag3 = t;
		int t_l = a1size;
		a1size = a2size;
		a2size = a3size;
		a3size = t_l;
		offset1 = offset2;
		offset2 = offset3;
		offset3 = minCol-1;
		__shared__ short temp[N_THREADS];
		initAntiDiag3(antiDiag3, a3size, offset3, maxCol, antiDiagNo, best - scoreDropOff, GAP_EXT, UNDEF);
		
		computeAntidiag(antiDiag1, antiDiag2, antiDiag3, querySeg, databaseSeg, best, scoreDropOff, cols, rows, minCol, maxCol, antiDiagNo, offset1, offset2, direction);	 	
		__syncthreads();	
	
		int tmp, antiDiagBest = UNDEF;	
		for(int i=0; i<a3size; i+=N_THREADS){
			int size = a3size-i;
			
			if(myTId<N_THREADS){
				temp[myTId] = (myTId<size) ? antiDiag3[myTId+i]:UNDEF;				
			}
			__syncthreads();
			
			tmp = reduce_max(temp,size);
			antiDiagBest = (tmp>antiDiagBest) ? tmp:antiDiagBest;

		}
		best = (best > antiDiagBest) ? best : antiDiagBest;
		//int prova = simple_max(antiDiag3, a3size);	
		//if(prova!=antiDiagBest){
		//	if(myTId==0)
		//		printf("errore %d/%d\n", prova,antiDiagBest);
		//}
		while (minCol - offset3 < a3size && antiDiag3[minCol - offset3] == UNDEF &&
			   minCol - offset2 - 1 < a2size && antiDiag2[minCol - offset2 - 1] == UNDEF)
		{
			++minCol;
		}

		// Calculate new maxCol
		while (maxCol - offset3 > 0 && (antiDiag3[maxCol - offset3 - 1] == UNDEF) &&
									   (antiDiag2[maxCol - offset2 - 1] == UNDEF))
		{
			--maxCol;
		}
		++maxCol;

		// Calculate new lowerDiag and upperDiag of extended seed
		calcExtendedLowerDiag(lowerDiag, minCol, antiDiagNo);
		calcExtendedUpperDiag(upperDiag, maxCol - 1, antiDiagNo);
		
		// end of databaseSeg reached?
		minCol = (minCol > (antiDiagNo + 2 - rows)) ? minCol : (antiDiagNo + 2 - rows);
		// end of querySeg reached?
		maxCol = (maxCol < cols) ? maxCol : cols;
		//}
	}

	int longestExtensionCol = a3size + offset3 - 2;
	int longestExtensionRow = antiDiagNo - longestExtensionCol;
	int longestExtensionScore = antiDiag3[longestExtensionCol - offset3];
	
	if (longestExtensionScore == UNDEF)
	{
		if (antiDiag2[a2size -2] != UNDEF)
		{
			// reached end of query segment
			longestExtensionCol = a2size + offset2 - 2;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
			
		}
		else if (a2size > 2 && antiDiag2[a2size-3] != UNDEF)
		{
			// reached end of database segment
			longestExtensionCol = a2size + offset2 - 3;
			longestExtensionRow = antiDiagNo - 1 - longestExtensionCol;
			longestExtensionScore = antiDiag2[longestExtensionCol - offset2];
			
		}
	}


	//could be parallelized in some way
	if (longestExtensionScore == UNDEF){

		// general case
		for (int i = 0; i < a1size; ++i){

			if (antiDiag1[i] > longestExtensionScore){

				longestExtensionScore = antiDiag1[i];
				longestExtensionCol = i + offset1;
				longestExtensionRow = antiDiagNo - 2 - longestExtensionCol;
			
			}
		}
	}
	
	if (longestExtensionScore != UNDEF)
		updateExtendedSeedL(mySeed, direction, longestExtensionCol, longestExtensionRow, lowerDiag, upperDiag);
	seed[myId] = mySeed;
	res[myId] = longestExtensionScore;
	//}

}

inline void extendSeedL(vector<SeedL> &seeds,
			ExtensionDirectionL direction,
			vector<string> &target,
			vector<string> &query,
			vector<ScoringSchemeL> &penalties,
			int const& XDrop,
			int const& kmer_length,
			int *res,
			int numAlignments,
			int ngpus
			)
{

	if(scoreGapExtend(penalties[0]) >= 0){

		cout<<"Error: Logan does not support gap extension penalty >= 0\n";
		exit(-1);
	}
	if(scoreGapOpen(penalties[0]) >= 0){

		cout<<"Error: Logan does not support gap opening penalty >= 0\n";
		exit(-1);
	}
	//start measuring time
	auto start_t1 = NOW;

	//declare streams
	cudaStream_t stream_r[MAX_GPUS], stream_l[MAX_GPUS];

	// NB nSequences is correlated to the number of GPUs that we have
	
	int nSequences = numAlignments/ngpus;
	//int nSeqInt = nSequences*sizeof(int);
	int nSequencesLast = nSequences+numAlignments%ngpus;
	//int nSeqIntLast = nSequencesLast*sizeof(int);	

	//final result of the alignment
	int *scoreLeft = (int *)malloc(numAlignments * sizeof(int));
	int *scoreRight = (int *)malloc(numAlignments * sizeof(int));

	//create two sets of seeds
	//copy seeds
	vector<SeedL> seeds_r;
	vector<SeedL> seeds_l;
	seeds_r.reserve(numAlignments);
	//seeds_l.reserve(numAlignments);
	for (int i=0; i<seeds.size(); i++){
			seeds_r.push_back(seeds[i]);
			//seeds_l.push_back(seeds[i]);
	}

	//sequences offsets	 		
	vector<int> offsetLeftQ[MAX_GPUS];
	vector<int> offsetLeftT[MAX_GPUS];	
	vector<int> offsetRightQ[MAX_GPUS];	
	vector<int> offsetRightT[MAX_GPUS];

	//shared_mem_size per block per GPU
	int shared_left[MAX_GPUS];
	int shared_right[MAX_GPUS];

	//antidiag in case shared memory isn't enough
	short *ant_l[MAX_GPUS], *ant_r[MAX_GPUS];

	//boolean used to know which kernel will be launched
	bool global_left[MAX_GPUS], global_right[MAX_GPUS];

	//total lenght of the sequences
	int totalLengthQPref[MAX_GPUS];
	int totalLengthTPref[MAX_GPUS];
	int totalLengthQSuff[MAX_GPUS];
	int totalLengthTSuff[MAX_GPUS];

	//declare and allocate sequences prefixes and suffixes
	char *prefQ[MAX_GPUS], *prefT[MAX_GPUS];
	char *suffQ[MAX_GPUS], *suffT[MAX_GPUS];

	//declare GPU offsets
	int *offsetLeftQ_d[MAX_GPUS], *offsetLeftT_d[MAX_GPUS];
	int *offsetRightQ_d[MAX_GPUS], *offsetRightT_d[MAX_GPUS];
	
	//declare GPU results
	int *scoreLeft_d[MAX_GPUS], *scoreRight_d[MAX_GPUS];

	//declare GPU seeds
	SeedL *seed_d_l[MAX_GPUS], *seed_d_r[MAX_GPUS];

	//declare prefixes and suffixes on the GPU  
	char *prefQ_d[MAX_GPUS], *prefT_d[MAX_GPUS];
	char *suffQ_d[MAX_GPUS], *suffT_d[MAX_GPUS];

	#pragma omp parallel for
	for(int i = 0; i < ngpus; i++){
		int dim = nSequences;
		if(i==ngpus-1)
			dim = nSequencesLast;
		//compute offsets and shared memory per block
		shared_left[i]=0;
		shared_right[i]=0;
		for(int j = 0; j < dim; j++){

			offsetLeftQ[i].push_back(getBeginPositionV(seeds[j+i*nSequences]));
			offsetLeftT[i].push_back(getBeginPositionH(seeds[j+i*nSequences]));
			shared_left[i] = std::max(std::min(offsetLeftQ[i][j],offsetLeftT[i][j]), shared_left[i]);
			
			offsetRightQ[i].push_back(query[j+i*nSequences].size()-getEndPositionV(seeds[j+i*nSequences]));
			offsetRightT[i].push_back(target[j+i*nSequences].size()-getEndPositionH(seeds[j+i*nSequences]));
			shared_right[i] = std::max(std::min(offsetRightQ[i][j], offsetRightT[i][j]), shared_right[i]);
		}
		//check if shared memory is enough
		global_left[i] = false;
		global_right[i] = false;
		if(shared_left[i]>=MAX_SIZE_ANTIDIAG){
			cudaErrchk(cudaMalloc(&ant_l[i], sizeof(short)*shared_left[i]*3*dim));
			global_left[i] = true;
			cout<<"LEFT GLOBAL GPU "<< i <<endl;
		}
		if(shared_right[i]>=MAX_SIZE_ANTIDIAG){
			cudaErrchk(cudaMalloc(&ant_r[i], sizeof(short)*shared_right[i]*3*dim));
			global_right[i] = true;
			cout<<"RIGHT GLOBAL GPU "<< i <<endl;
		}
		//compute antidiagonal offsets
		partial_sum(offsetLeftQ[i].begin(),offsetLeftQ[i].end(),offsetLeftQ[i].begin());	
		partial_sum(offsetLeftT[i].begin(),offsetLeftT[i].end(),offsetLeftT[i].begin());
		partial_sum(offsetRightQ[i].begin(),offsetRightQ[i].end(),offsetRightQ[i].begin());
		partial_sum(offsetRightT[i].begin(),offsetRightT[i].end(),offsetRightT[i].begin());
		//set total length of the sequences
		totalLengthQPref[i] = offsetLeftQ[i][dim-1];
		totalLengthTPref[i] = offsetLeftT[i][dim-1];
		totalLengthQSuff[i] = offsetRightQ[i][dim-1];
		totalLengthTSuff[i] = offsetRightT[i][dim-1];
		//allocate sequences prefix and suffix on the CPU
		prefQ[i] = (char*)malloc(sizeof(char)*totalLengthQPref[i]);
		prefT[i] = (char*)malloc(sizeof(char)*totalLengthTPref[i]);
		suffQ[i] = (char*)malloc(sizeof(char)*totalLengthQSuff[i]);
		suffT[i] = (char*)malloc(sizeof(char)*totalLengthTSuff[i]);
		//generate prefix and suffix on the CPU
		//std::cout << "SETTING UP PREF/SUFF" << std::endl;
		reverse_copy(query[0+i*nSequences].c_str(),query[0+i*nSequences].c_str()+offsetLeftQ[i][0],prefQ[i]);
		//memcpy(prefQ[i], query[0+i*nSequences].c_str(), offsetLeftQ[i][0]);
		memcpy(prefT[i], target[0+i*nSequences].c_str(), offsetLeftT[i][0]);
		memcpy(suffQ[i], query[0+i*nSequences].c_str()+getEndPositionV(seeds[0+i*nSequences]), offsetRightQ[i][0]);
		reverse_copy(target[0+i*nSequences].c_str()+getEndPositionH(seeds[0+i*nSequences]),target[0+i*nSequences].c_str()+getEndPositionH(seeds[0+i*nSequences])+offsetRightT[i][0],suffT[i]);
		//memcpy(suffT[i], target[0+i*nSequences].c_str()+getEndPositionH(seeds[0+i*nSequences]), offsetRightT[i][0]);
		for(int j = 1; j<dim; j++){
			char *seqptr = prefQ[i] + offsetLeftQ[i][j-1];
			reverse_copy(query[j+i*nSequences].c_str(),query[j+i*nSequences].c_str()+(offsetLeftQ[i][j]-offsetLeftQ[i][j-1]),seqptr);
			//memcpy(seqptr, query[j+i*nSequences].c_str(), offsetLeftQ[i][j]-offsetLeftQ[i][j-1]);
			seqptr = prefT[i] + offsetLeftT[i][j-1];
			memcpy(seqptr, target[j+i*nSequences].c_str(), offsetLeftT[i][j]-offsetLeftT[i][j-1]);
			seqptr = suffQ[i] + offsetRightQ[i][j-1];
			memcpy(seqptr, query[j+i*nSequences].c_str()+getEndPositionV(seeds[j+i*nSequences]), offsetRightQ[i][j]-offsetRightQ[i][j-1]);
			seqptr = suffT[i] + offsetRightT[i][j-1];
			reverse_copy(target[j+i*nSequences].c_str()+getEndPositionH(seeds[j+i*nSequences]),target[j+i*nSequences].c_str()+getEndPositionH(seeds[j+i*nSequences])+(offsetRightT[i][j]-offsetRightT[i][j-1]),seqptr);
			//memcpy(seqptr, target[j+i*nSequences].c_str()+getEndPositionH(seeds[j+i*nSequences]), offsetRightT[i][j]-offsetRightT[i][j-1]);

		}
	}

	for(int i = 0; i < ngpus; i++){
		int dim = nSequences;
		if(i==ngpus-1)
			dim = nSequencesLast;
		//set gpu device
		cudaSetDevice(i);
		//create streams
		cudaStreamCreateWithFlags(&stream_r[i],cudaStreamNonBlocking);
		cudaStreamCreateWithFlags(&stream_l[i],cudaStreamNonBlocking);
		//allocate offsets on the GPU
		cudaErrchk(cudaMalloc(&offsetLeftQ_d[i], dim*sizeof(int)));
		cudaErrchk(cudaMalloc(&offsetLeftT_d[i], dim*sizeof(int)));
		cudaErrchk(cudaMalloc(&offsetRightQ_d[i], dim*sizeof(int)));
		cudaErrchk(cudaMalloc(&offsetRightT_d[i], dim*sizeof(int)));
		//allocate results on the GPU
		cudaErrchk(cudaMalloc(&scoreLeft_d[i], dim*sizeof(int)));
		cudaErrchk(cudaMalloc(&scoreRight_d[i], dim*sizeof(int)));
		//allocate seeds on the GPU
		cudaErrchk(cudaMalloc(&seed_d_l[i], dim*sizeof(SeedL)));
		cudaErrchk(cudaMalloc(&seed_d_r[i], dim*sizeof(SeedL)));
		//allocate sequences on the GPU
		cudaErrchk(cudaMalloc(&prefQ_d[i], totalLengthQPref[i]*sizeof(char)));
		cudaErrchk(cudaMalloc(&prefT_d[i], totalLengthTPref[i]*sizeof(char)));
		cudaErrchk(cudaMalloc(&suffQ_d[i], totalLengthQSuff[i]*sizeof(char)));
		cudaErrchk(cudaMalloc(&suffT_d[i], totalLengthTSuff[i]*sizeof(char)));
		//copy seeds on the GPU
		cudaErrchk(cudaMemcpyAsync(seed_d_l[i], &seeds[0]+i*nSequences, dim*sizeof(SeedL), cudaMemcpyHostToDevice, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(seed_d_r[i], &seeds_r[0]+i*nSequences, dim*sizeof(SeedL), cudaMemcpyHostToDevice, stream_r[i]));
		//copy offsets on the GPU
		cudaErrchk(cudaMemcpyAsync(offsetLeftQ_d[i], &offsetLeftQ[i][0], dim*sizeof(int), cudaMemcpyHostToDevice, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(offsetLeftT_d[i], &offsetLeftT[i][0], dim*sizeof(int), cudaMemcpyHostToDevice, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(offsetRightQ_d[i], &offsetRightQ[i][0], dim*sizeof(int), cudaMemcpyHostToDevice, stream_r[i]));
		cudaErrchk(cudaMemcpyAsync(offsetRightT_d[i], &offsetRightT[i][0], dim*sizeof(int), cudaMemcpyHostToDevice, stream_r[i]));
		//copy sequences on the GPU
		cudaErrchk(cudaMemcpyAsync(prefQ_d[i], prefQ[i], totalLengthQPref[i]*sizeof(char), cudaMemcpyHostToDevice, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(prefT_d[i], prefT[i], totalLengthTPref[i]*sizeof(char), cudaMemcpyHostToDevice, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(suffQ_d[i], suffQ[i], totalLengthQSuff[i]*sizeof(char), cudaMemcpyHostToDevice, stream_r[i]));
		cudaErrchk(cudaMemcpyAsync(suffT_d[i], suffT[i], totalLengthTSuff[i]*sizeof(char), cudaMemcpyHostToDevice, stream_r[i]));
		//OK
	}
	
	
	auto end_t1 = NOW;
	duration<double> transfer1=end_t1-start_t1;
	std::cout << "Input setup time: " << transfer1.count() << std::endl;
	auto start_c = NOW;
	
	//execute kernels
	for(int i = 0; i<ngpus;i++){
		cudaSetDevice(i);
		
		int dim = nSequences;
		if(i==ngpus-1)
			dim = nSequencesLast;
		
		if(global_left[i]){
			extendSeedLGappedXDropOneDirectionGlobal <<<dim, N_THREADS, 0, stream_l[i]>>> (seed_d_l[i], prefQ_d[i], prefT_d[i], EXTEND_LEFTL, XDrop, scoreLeft_d[i], offsetLeftQ_d[i], offsetLeftT_d[i], shared_left[i], ant_l[i]);
		}else{
			extendSeedLGappedXDropOneDirectionShared <<<dim, N_THREADS, 3*shared_left[i]*sizeof(short), stream_l[i]>>> (seed_d_l[i], prefQ_d[i], prefT_d[i], EXTEND_LEFTL, XDrop, scoreLeft_d[i], offsetLeftQ_d[i], offsetLeftT_d[i], shared_left[i]);
		}
		if(global_right[i]){
			extendSeedLGappedXDropOneDirectionGlobal <<<dim, N_THREADS, 0, stream_r[i]>>> (seed_d_r[i], suffQ_d[i], suffT_d[i], EXTEND_RIGHTL, XDrop, scoreRight_d[i], offsetRightQ_d[i], offsetRightT_d[i], shared_right[i], ant_r[i]);
		}else{ 
			extendSeedLGappedXDropOneDirectionShared <<<dim, N_THREADS, 3*shared_right[i]*sizeof(short), stream_r[i]>>> (seed_d_r[i], suffQ_d[i], suffT_d[i], EXTEND_RIGHTL, XDrop, scoreRight_d[i], offsetRightQ_d[i], offsetRightT_d[i], shared_right[i]);
		}

		//cout<<"LAUNCHED"<<endl;
	}
	
	for(int i = 0; i < ngpus; i++){
		cudaSetDevice(i);
		int dim = nSequences;
		if(i==ngpus-1)
			dim = nSequencesLast;
		cudaErrchk(cudaMemcpyAsync(scoreLeft+i*nSequences, scoreLeft_d[i], dim*sizeof(int), cudaMemcpyDeviceToHost, stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(&seeds[0]+i*nSequences, seed_d_l[i], dim*sizeof(SeedL), cudaMemcpyDeviceToHost,stream_l[i]));
		cudaErrchk(cudaMemcpyAsync(scoreRight+i*nSequences, scoreRight_d[i], dim*sizeof(int), cudaMemcpyDeviceToHost, stream_r[i]));
		cudaErrchk(cudaMemcpyAsync(&seeds_r[0]+i*nSequences, seed_d_r[i], dim*sizeof(SeedL), cudaMemcpyDeviceToHost,stream_r[i]));
	}

	for(int i = 0; i < ngpus; i++){
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}

	auto end_c = NOW;
	duration<double> compute = end_c-start_c;
	std::cout << "Compute time: " << compute.count() << std::endl;

	cudaErrchk(cudaPeekAtLastError());

	//cudaStreamDestroy(stream_l);
	//cudaStreamDestroy(stream_r);
	auto start_f = NOW;


	for(int i = 0; i < ngpus; i++){
		cudaSetDevice(i);

		cudaStreamDestroy(stream_l[i]);
		cudaStreamDestroy(stream_r[i]);
		free(prefQ[i]);
		free(prefT[i]);
		free(suffQ[i]);
		free(suffT[i]);
		cudaErrchk(cudaFree(prefQ_d[i]));
		cudaErrchk(cudaFree(prefT_d[i]));
		cudaErrchk(cudaFree(suffQ_d[i]));
		cudaErrchk(cudaFree(suffT_d[i]));
		cudaErrchk(cudaFree(offsetLeftQ_d[i]));
		cudaErrchk(cudaFree(offsetLeftT_d[i]));
		cudaErrchk(cudaFree(offsetRightQ_d[i]));
		cudaErrchk(cudaFree(offsetRightT_d[i]));
		cudaErrchk(cudaFree(seed_d_l[i]));
		cudaErrchk(cudaFree(seed_d_r[i]));
		cudaErrchk(cudaFree(scoreLeft_d[i]));
		cudaErrchk(cudaFree(scoreRight_d[i]));
	
		if(global_left[i])
			cudaErrchk(cudaFree(ant_l[i])); 
			if(global_right[i])
			cudaErrchk(cudaFree(ant_r[i]));

	}
	auto end_f = NOW;
	
	for(int i = 0; i < numAlignments; i++){
		res[i] = scoreLeft[i]+scoreRight[i]+kmer_length;
		setEndPositionH(seeds[i], getEndPositionH(seeds_r[i]));    
		setEndPositionV(seeds[i], getEndPositionV(seeds_r[i])); 
		//cout << res[i] <<endl;
	}
	
	free(scoreLeft);
	free(scoreRight);
		
}

