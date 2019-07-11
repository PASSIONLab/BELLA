//==================================================================
// Title:  C++ x-drop seed-and-extend alignment algorithm
// Author: G. Guidi, A. Zeni
// Date:   6 March 2019
//==================================================================

#include<algorithm> 
#include<cassert>

template<typename Tx_>
const Tx_&  min(const Tx_& _Left, const Tx_& Right_)
{   // return smaller of _Left and Right_
    if (_Left < Right_)
        return _Left;
    else
        return Right_;
}

template<typename Tx_, typename Ty_>
Tx_  min(const Tx_& _Left, const Ty_& Right_)
{   // return smaller of _Left and Right_
    return (Right_ < _Left ? Right_ : _Left);
}
template<typename Ty_>
Ty_ const &
__device__ __host__ max_logan(const Ty_& _Left, const Ty_& Right_)
{   // return larger of _Left and Right_
    if (_Left < Right_)
        return Right_;
    else
        return _Left;
}

template<typename Tx_, typename Ty_>
Tx_
__device__ __host__ max_logan(const Tx_& _Left, const Ty_& Right_)
{   // return smaller of _Left and Right_
    return (Right_ < _Left ? _Left : Right_);
}


struct SeedL
{
	int beginPositionH;
	int beginPositionV;
	int endPositionH;
	int endPositionV;
	int seedLength;
	int lowerDiagonal;  // GGGG: it might possibly be a std::string
	int upperDiagonal;  // GGGG: it might possibly be a std::string
	int beginDiagonal;
	int endDiagonal;
	int score;

	SeedL(): beginPositionH(0), beginPositionV(0), endPositionH(0), endPositionV(0), lowerDiagonal(0), upperDiagonal(0), score(0)
	{}

	SeedL(int beginPositionH, int beginPositionV, int seedLength):
		beginPositionH(beginPositionH), beginPositionV(beginPositionV), endPositionH(beginPositionH + seedLength),
		endPositionV(beginPositionV + seedLength), lowerDiagonal((beginPositionH - beginPositionV)),
		upperDiagonal((beginPositionH - beginPositionV)), beginDiagonal(beginPositionH - beginPositionV),
		endDiagonal(endPositionH - endPositionV), score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedL(int beginPositionH, int beginPositionV, int endPositionH, int endPositionV):
		beginPositionH(beginPositionH),
		beginPositionV(beginPositionV),
		endPositionH(endPositionH),
		endPositionV(endPositionV),
		lowerDiagonal(min((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		upperDiagonal(max((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		beginDiagonal((beginPositionH - beginPositionV)),
		endDiagonal((endPositionH - endPositionV)),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	__device__ __host__ SeedL(SeedL const& other):
		beginPositionH(other.beginPositionH),
		beginPositionV(other.beginPositionV),
		endPositionH(other.endPositionH),
		endPositionV(other.endPositionV),
		lowerDiagonal(other.lowerDiagonal),
		upperDiagonal(other.upperDiagonal),
		beginDiagonal(other.beginDiagonal),
		endDiagonal(other.endDiagonal),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

};

struct Result
{
	SeedL myseed;
	int score; 			// alignment score
	int length;			// overlap length / max extension

	Result() : score(0), length(0)//check
	{
		//myseed=SeedL();
	}

	Result(int kmerLen) : score(0), length(kmerLen)
	{
		//myseed=SeedL();
	}

};

// GGGG we can think about this later
// AAAA add setter also

int
__device__ getAlignScore(SeedL const &myseed){
	return myseed.score;
}

int
__device__ __host__ getBeginPositionH(SeedL const &myseed){
	return myseed.beginPositionH;
}

int
__device__ __host__ getBeginPositionV(SeedL const &myseed){
	return myseed.beginPositionV;
}

int
__device__ __host__ getEndPositionH(SeedL const &myseed){
	return myseed.endPositionH;
}

int
__device__ __host__ getEndPositionV(SeedL const &myseed){
	return myseed.endPositionV;
}

int
__device__ getSeedLLength(SeedL const &myseed){
	return myseed.seedLength;
}

int
__device__ getLowerDiagonal(SeedL const &myseed){
	return myseed.lowerDiagonal;
}

int
__device__ getUpperDiagonal(SeedL const &myseed){
	return myseed.upperDiagonal;
}

int
__device__ getBeginDiagonal(SeedL const &myseed){
	return myseed.beginDiagonal;
}

int
__device__ getEndDiagonal(SeedL const &myseed){
	return myseed.endDiagonal;
}

void
__device__ setAlignScore(SeedL &myseed,int const value){
	myseed.score = value;
}

void
__device__ __host__ setBeginPositionH(SeedL &myseed,int const value){
	myseed.beginPositionH = value;
}

void
__device__ __host__ setBeginPositionV(SeedL &myseed,int const value){
	myseed.beginPositionV = value;
}

void
__device__ __host__ setEndPositionH(SeedL &myseed,int const value){
	myseed.endPositionH = value;
}

void
__device__ __host__ setEndPositionV(SeedL &myseed,int const value){
	myseed.endPositionV = value;
}

void
__device__ __host__ setSeedLLength(SeedL &myseed,int const value){
	myseed.seedLength = value;
}

void
__device__ __host__ setLowerDiagonal(SeedL &myseed,int const value){
	myseed.lowerDiagonal = value;
}

void
__device__ __host__ setUpperDiagonal(SeedL &myseed,int const value){
	myseed.upperDiagonal = value;
}

void
__device__ __host__ setBeginDiagonal(SeedL &myseed,int const value){
	myseed.beginDiagonal = value;
}

void
__device__ __host__ setEndDiagonal(SeedL &myseed,int const value){
	myseed.endDiagonal = value;
}
