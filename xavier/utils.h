//===========================================================================
// Title:  Xavier: High-Performance X-Drop Adaptive Banded Pairwise Alignment
// Author: G. Guidi, E. Younis, A. Zeni
// Date:   8 March 2019
//===========================================================================

#ifndef __NVCC__

#ifndef UTILS_H
#define UTILS_H

#include<algorithm> 
#include<cassert>

//===========================================================================
// UTILS
//===========================================================================

struct SeedX
{
	int beginPositionH;
	int beginPositionV;
	int endPositionH;
	int endPositionV;
	int SeedXength;
	int lowerDiagonal;
	int upperDiagonal;
	int beginDiagonal;
	int endDiagonal;
	int score;

	SeedX(): beginPositionH(0), beginPositionV(0), endPositionH(0), endPositionV(0), lowerDiagonal(0), upperDiagonal(0), score(0)
	{}

	SeedX(int beginPositionH, int beginPositionV, int SeedXength):
		beginPositionH(beginPositionH), beginPositionV(beginPositionV), endPositionH(beginPositionH + SeedXength),
		endPositionV(beginPositionV + SeedXength), lowerDiagonal((beginPositionH - beginPositionV)),
		upperDiagonal((beginPositionH - beginPositionV)), beginDiagonal(beginPositionH - beginPositionV),
		endDiagonal(endPositionH - endPositionV), score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedX(int beginPositionH, int beginPositionV, int endPositionH, int endPositionV):
		beginPositionH(beginPositionH),
		beginPositionV(beginPositionV),
		endPositionH(endPositionH),
		endPositionV(endPositionV),
		lowerDiagonal(std::min((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		upperDiagonal(std::max((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		beginDiagonal((beginPositionH - beginPositionV)),
		endDiagonal((endPositionH - endPositionV)),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedX(SeedX const& other):
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

inline int
getAlignScore(SeedX const &myseed){
	return myseed.score;
}

inline int
getBeginPositionH(SeedX const &myseed){
	return myseed.beginPositionH;
}

inline int
getBeginPositionV(SeedX const &myseed){
	return myseed.beginPositionV;
}

inline int
getEndPositionH(SeedX const &myseed){
	return myseed.endPositionH;
}

inline int
getEndPositionV(SeedX const &myseed){
	return myseed.endPositionV;
}

inline int
getSeedXLength(SeedX const &myseed){
	return myseed.SeedXength;
}

inline int
getLowerDiagonal(SeedX const &myseed){
	return myseed.lowerDiagonal;
}

inline int
getUpperDiagonal(SeedX const &myseed){
	return myseed.upperDiagonal;
}

inline int
getBeginDiagonal(SeedX const &myseed){
	return myseed.beginDiagonal;
}

inline int
getEndDiagonal(SeedX const &myseed){
	return myseed.endDiagonal;
}

inline void
setAlignScore(SeedX &myseed,int const value){
	myseed.score = value;
}

inline void
setBeginPositionH(SeedX &myseed,int const value){
	myseed.beginPositionH = value;
}

inline void
setBeginPositionV(SeedX &myseed,int const value){
	myseed.beginPositionV = value;
}

inline void
setEndPositionH(SeedX &myseed,int const value){
	myseed.endPositionH = value;
}

inline void
setEndPositionV(SeedX &myseed,int const value){
	myseed.endPositionV = value;
}

inline void
setSeedXLength(SeedX &myseed,int const value){
	myseed.SeedXength = value;
}

inline void
setLowerDiagonal(SeedX &myseed,int const value){
	myseed.lowerDiagonal = value;
}

inline void
setUpperDiagonal(SeedX &myseed,int const value){
	myseed.upperDiagonal = value;
}

inline void
setBeginDiagonal(SeedX &myseed,int const value){
	myseed.beginDiagonal = value;
}

inline void
setEndDiagonal(SeedX &myseed,int const value){
	myseed.endDiagonal = value;
}

#endif

#endif // __NVCC__