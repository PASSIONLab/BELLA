//==================================================================
// Title:  LOGAN: X-Drop Adaptive Banded Alignment
// Author: G. Guidi, E. Younis
// Date:   30 April 2019
//==================================================================

#ifndef SIMD_UTILS_H
#define SIMD_UTILS_H

#include <cstdint>
#include"score.h"
#include"utils.h"

//======================================================================================
// GLOBAL FUNCTION DECLARATION
//======================================================================================

#ifndef __AVX2__
#define __AVX2__
#endif

#ifdef  __AVX2__ 	// Compile flag: -mavx2
#define VECTORWIDTH  (32)
#define LOGICALWIDTH (VECTORWIDTH - 1)
#define vectorType   __m256i
#define addOp    	_mm256_adds_epi8  	// saturated arithmetic
#define subOp    	_mm256_subs_epi8  	// saturated arithmetic
#define maxOp    	_mm256_max_epi8   	// max
#define setOp   	_mm256_set1_epi8  	// set1 operation
#define blendvOp	_mm256_blendv_epi8  // blending operation
#define cmpeqOp 	_mm256_cmpeq_epi8 	// compare equality operation
#elif __SSE4_2__ 	// Compile flag: -msse4.2
#define VECTORWIDTH  (8)
#define LOGICALWIDTH (VECTORWIDTH - 1)
#define vectorType  __m128i
#define addOp    	_mm_adds_epi16 	// saturated arithmetic
#define subOp    	_mm_subs_epi16  // saturated arithmetic
#define maxOp    	_mm_max_epi16  	// max
#define setOp   	_mm_set1_epi16  // set1 operation
#define blendvOp 	_mm_blendv_epi8 // blending operation
#define cmpeqOp  	_mm_cmpeq_epi16 // compare equality operation
#endif

//======================================================================================
// GLOBAL VARIABLE DEFINITION
//======================================================================================

#define NINF  	(std::numeric_limits<int8_t>::min())
#define goRIGHT (0)
#define goDOWN  (1)
#define MIDDLE 	(LOGICALWIDTH / 2)

#define CUTOFF	(std::numeric_limits<int8_t>::max() - 25)

#ifdef DEBUG
	#define myLog( var ) do { std::cerr << "LOG:	" << __FILE__ << "(" << __LINE__ << ")	" << #var << " = " << (var) << std::endl; } while(0)
#else
	#define myLog( var )
#endif

//======================================================================================
// SIMD UTILS
//======================================================================================

typedef int8_t element_t;

typedef union {
	vectorType  simd;
	element_t 	elem[VECTORWIDTH];
} vectorUnionType;

void
printVectorC(vectorType a) {

	vectorUnionType tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		printf("%c,", tmp.elem[i]);
	printf("%c}\n", tmp.elem[VECTORWIDTH-1]);
}

void
printVectorD(vectorType a) {

	vectorUnionType tmp;
	tmp.simd = a;

	printf("{");
	for (int i = 0; i < VECTORWIDTH-1; ++i)
		printf("%d,", tmp.elem[i]);
	printf("%d}\n", tmp.elem[VECTORWIDTH-1]);
}

enum ExtDirectionL
{
	LOGAN_EXTEND_NONE  = 0,
	LOGAN_EXTEND_LEFT  = 1,
	LOGAN_EXTEND_RIGHT = 2,
	LOGAN_EXTEND_BOTH  = 3
};

// TODO: REWRITE shiftLeft/RightShift for int8_t
// TODO: Verify
#ifdef __AVX2__
inline vectorUnionType
shiftLeft (const vectorType& _a) { // this work for avx2

	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 1);
	b.elem[VECTORWIDTH - 1] = NINF;

	return b;
}
inline vectorUnionType
shiftRight (const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 1);
	b.elem[0] = NINF;
	return b;
}
#elif __SSE4_2__
inline vectorUnionType
shiftLeft(const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(_mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(2, 0, 0, 1)), a.simd, 2);
	b.elem[VECTORWIDTH - 1] = NINF;
	return b;
}
inline vectorUnionType
shiftRight(const vectorType& _a) { // this work for avx2
	vectorUnionType a;
	a.simd = _a;

	vectorUnionType b;
	// https://stackoverflow.com/questions/25248766/emulating-shifts-on-32-bytes-with-avx
	b.simd = _mm256_alignr_epi8(a.simd, _mm256_permute2x128_si256(a.simd, a.simd, _MM_SHUFFLE(0, 0, 2, 0)), 16 - 2);
	b.elem[0] = NINF;
	return b;
}
#endif

class LoganState
{
public:
	LoganState
	(
		SeedL& _seed,
	 	std::string const& targetSeg,
		std::string const& querySeg,
		ScoringSchemeL& scoringScheme,
		int64_t const &_scoreDropOff
	)
	{
		seed = _seed;

		hlength = targetSeg.length() + 1;
		vlength = querySeg.length()  + 1;

		if (hlength <= 1 || vlength <= 1)
		{
			std::string ErrorMessage = "Read(s) has length smaller than 1.";
			myLog(ErrorMessage);
			return;
		}

		// Convert from string to int array
		// This is the entire sequences
		queryh = new int8_t[hlength];
		queryv = new int8_t[vlength];
		std::copy(targetSeg.begin(), targetSeg.end(), queryh);
		std::copy(querySeg.begin(), querySeg.end(), queryv);

		// pay attention
		// myLog( "Verify Here (ASSERT)" );

		matchCost    = scoreMatch(scoringScheme   );
		mismatchCost = scoreMismatch(scoringScheme);
		gapCost      = scoreGap(scoringScheme     );

		vmatchCost    = setOp (matchCost   );
		vmismatchCost = setOp (mismatchCost);
		vgapCost      = setOp (gapCost     );
		vzeros        = _mm256_setzero_si256();

		hoffset = LOGICALWIDTH;
		voffset = LOGICALWIDTH;

		bestScore    = 0;
		currScore    = 0;
		scoreOffset  = 0;
		scoreDropOff = _scoreDropOff;
		xDropCond   = false;
	}

	~LoganState()
	{
		delete [] queryh;
		delete [] queryv;
	}

	// i think this can be smaller than 64bit
	int64_t get_score_offset  ( void ) { return scoreOffset;  }
	int64_t get_best_score    ( void ) { return bestScore;    }
	int64_t get_curr_score    ( void ) { return currScore;    }
	int64_t get_score_dropoff ( void ) { return scoreDropOff; }

	void set_score_offset ( int64_t _scoreOffset ) { scoreOffset = _scoreOffset; }
	void set_best_score   ( int64_t _bestScore   ) { bestScore   = _bestScore;   }
	void set_curr_score   ( int64_t _currScore   ) { currScore   = _currScore;   }

	int8_t get_match_cost    ( void ) { return matchCost;    }
	int8_t get_mismatch_cost ( void ) { return mismatchCost; }
	int8_t get_gap_cost      ( void ) { return gapCost;      }

	vectorType get_vqueryh ( void ) { return vqueryh.simd; }
	vectorType get_vqueryv ( void ) { return vqueryv.simd; }

	vectorType get_antiDiag1 ( void ) { return antiDiag1.simd; }
	vectorType get_antiDiag2 ( void ) { return antiDiag2.simd; }
	vectorType get_antiDiag3 ( void ) { return antiDiag3.simd; }

	vectorType get_vmatchCost    ( void ) { return vmatchCost;    }
	vectorType get_vmismatchCost ( void ) { return vmismatchCost; }
	vectorType get_vgapCost      ( void ) { return vgapCost;      }
	vectorType get_vzeros        ( void ) { return vzeros;        }

	void update_vqueryh ( uint8_t idx, int8_t value ) { vqueryh.elem[idx] = value; }
	void update_vqueryv ( uint8_t idx, int8_t value ) { vqueryv.elem[idx] = value; }

	void update_antiDiag1 ( uint8_t idx, int8_t value ) { antiDiag1.elem[idx] = value; }
	void update_antiDiag2 ( uint8_t idx, int8_t value ) { antiDiag2.elem[idx] = value; }
	void update_antiDiag3 ( uint8_t idx, int8_t value ) { antiDiag3.elem[idx] = value; }

	void broadcast_antiDiag1 ( int8_t value ) { antiDiag1.simd = setOp( value ); }
	void broadcast_antiDiag2 ( int8_t value ) { antiDiag2.simd = setOp( value ); }
	void broadcast_antiDiag3 ( int8_t value ) { antiDiag3.simd = setOp( value ); }

	void set_antiDiag1 ( vectorType vector ) { antiDiag1.simd = vector; }
	void set_antiDiag2 ( vectorType vector ) { antiDiag2.simd = vector; }
	void set_antiDiag3 ( vectorType vector ) { antiDiag3.simd = vector; }

	void moveRight ( void )
	{
		// (a) shift to the left on query horizontal
		vqueryh = shiftLeft( vqueryh.simd );
		vqueryh.elem[LOGICALWIDTH - 1] = queryh[hoffset++];

		// (b) shift left on updated vector 1
		// this places the right-aligned vector 2 as a left-aligned vector 1
		antiDiag1.simd = antiDiag2.simd;
		antiDiag1 = shiftLeft( antiDiag1.simd );
		antiDiag2.simd = antiDiag3.simd;
	}

	void moveDown ( void )
	{
		// (a) shift to the right on query vertical
		vqueryv = shiftRight( vqueryv.simd );
		vqueryv.elem[0] = queryv[voffset++];

		// (b) shift to the right on updated vector 2
		// this places the left-aligned vector 3 as a right-aligned vector 2
		antiDiag1.simd = antiDiag2.simd;
		antiDiag2.simd = antiDiag3.simd;
		antiDiag2 = shiftRight( antiDiag2.simd );
	}

	// Seed position (define starting position and need to be updated when exiting)
	SeedL seed;

	// Sequence Lengths
	unsigned int hlength;
	unsigned int vlength;

	// Sequences as ints
	int8_t* queryh;
	int8_t* queryv;

	// Sequence pointers
	int hoffset;
	int voffset;

	// Constant Scoring Values
	int8_t matchCost;
	int8_t mismatchCost;
	int8_t gapCost;

	// Constant Scoring Vectors
	vectorType vmatchCost;
	vectorType vmismatchCost;
	vectorType vgapCost;
	vectorType vzeros;

	// Computation Vectors
	vectorUnionType antiDiag1;
	vectorUnionType antiDiag2;
	vectorUnionType antiDiag3;

	vectorUnionType vqueryh;
	vectorUnionType vqueryv;

	// X-Drop Variables
	int64_t bestScore;
	int64_t currScore;
	int64_t scoreOffset;
	int64_t scoreDropOff;
	bool xDropCond;
};

void operator+=(LoganState& state1, const LoganState& state2)
{
	state1.bestScore = state1.bestScore + state2.bestScore;
	state1.currScore = state1.currScore + state2.currScore;
}

#endif