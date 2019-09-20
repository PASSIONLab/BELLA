//==================================================================
// Title:  C++ x-drop seed-and-extend alignment algorithm
// Author: G. Guidi, A. Zeni
// Date:   6 April 2019
//==================================================================

#ifndef UTILS_H
#define UTILS_H

#include<algorithm> 
#include<cassert>

//======================================================================================
// UTILS
//======================================================================================

struct SeedL
{
	int beginPositionH;
	int beginPositionV;
	int endPositionH;
	int endPositionV;
	int seedLength;
	int lowerDiagonal;
	int upperDiagonal;
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
		lowerDiagonal(std::min((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		upperDiagonal(std::max((beginPositionH - beginPositionV), (endPositionH - endPositionV))),
		beginDiagonal((beginPositionH - beginPositionV)),
		endDiagonal((endPositionH - endPositionV)),
		score(0)
	{
		assert(upperDiagonal >= lowerDiagonal);
	}

	SeedL(SeedL const& other):
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
getAlignScore(SeedL const &myseed){
	return myseed.score;
}

inline int
getBeginPositionH(SeedL const &myseed){
	return myseed.beginPositionH;
}

inline int
getBeginPositionV(SeedL const &myseed){
	return myseed.beginPositionV;
}

inline int
getEndPositionH(SeedL const &myseed){
	return myseed.endPositionH;
}

inline int
getEndPositionV(SeedL const &myseed){
	return myseed.endPositionV;
}

inline int
getSeedLLength(SeedL const &myseed){
	return myseed.seedLength;
}

inline int
getLowerDiagonal(SeedL const &myseed){
	return myseed.lowerDiagonal;
}

inline int
getUpperDiagonal(SeedL const &myseed){
	return myseed.upperDiagonal;
}

inline int
getBeginDiagonal(SeedL const &myseed){
	return myseed.beginDiagonal;
}

inline int
getEndDiagonal(SeedL const &myseed){
	return myseed.endDiagonal;
}

inline void
setAlignScore(SeedL &myseed,int const value){
	myseed.score = value;
}

inline void
setBeginPositionH(SeedL &myseed,int const value){
	myseed.beginPositionH = value;
}

inline void
setBeginPositionV(SeedL &myseed,int const value){
	myseed.beginPositionV = value;
}

inline void
setEndPositionH(SeedL &myseed,int const value){
	myseed.endPositionH = value;
}

inline void
setEndPositionV(SeedL &myseed,int const value){
	myseed.endPositionV = value;
}

inline void
setSeedLLength(SeedL &myseed,int const value){
	myseed.seedLength = value;
}

inline void
setLowerDiagonal(SeedL &myseed,int const value){
	myseed.lowerDiagonal = value;
}

inline void
setUpperDiagonal(SeedL &myseed,int const value){
	myseed.upperDiagonal = value;
}

inline void
setBeginDiagonal(SeedL &myseed,int const value){
	myseed.beginDiagonal = value;
}

inline void
setEndDiagonal(SeedL &myseed,int const value){
	myseed.endDiagonal = value;
}

#endif