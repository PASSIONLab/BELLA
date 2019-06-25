#ifndef BFG_KMER_HPP
#define BFG_KMER_HPP

#include <stdio.h>
#include <stdint.h>
#include <cassert>
#include <cstring>
#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <functional>
#include <cstdint>

#include "hash_funcs.h"
#include "common.h"



/* Short description: 
 *  - Store kmer strings by using 2 bits per base instead of 8 
 *  - Easily return reverse complements of kmers, e.g. TTGG -> CCAA
 *  - Easily compare kmers
 *  - Provide hash of kmers
 *  - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
 *  */
#define N_LONGS (MAX_KMER_SIZE/32)
#define N_BYTES (MAX_KMER_SIZE/4)

class Kmer {
 public:
  typedef std::array<uint64_t, N_LONGS> MERARR;
  typedef std::array<uint8_t, N_BYTES> BYTEARR;

  Kmer();
  Kmer(const Kmer& o);
  explicit Kmer(const char *s, unsigned int len);
  void copyDataFrom(uint8_t *  mybytes)	// this is like a shadow constructor (to avoid accidental signature match with the existing constructor)
  {
	  memcpy(longs.data(), mybytes, sizeof(uint64_t) * (N_LONGS));
  }
  explicit Kmer(const MERARR & arr)
  {
	std::memcpy (longs.data(), arr.data(), sizeof(uint64_t) * (N_LONGS));
  }
  static std::vector<Kmer> getKmers(std::string seq);

  Kmer& operator=(const Kmer& o);  
  void set_deleted();
  bool operator<(const Kmer& o) const;
  bool operator==(const Kmer& o) const;
  bool operator!=(const Kmer& o) const {
    return !(*this == o);
  }

  void set_kmer(const char *s, unsigned int len);
  uint64_t hash() const;
  

  Kmer twin() const;
  Kmer rep() const; // ABAB: return the smaller of itself (lexicographically) or its reversed-complement (i.e. twin)
  
  // NOTE: commented out for version with Kmers of variable lengths because of lack of testing
  // Kmer getLink(const size_t index) const;
  // Kmer forwardBase(const char b) const;
  // Kmer backwardBase(const char b) const;
  // bool equalUpToLastBase(const Kmer & rhs);	// returns true for completely identical k-mers as well as k-mers that only differ at the last base

  Kmer hopc() const;
  std::string getBinary() const;  
  void toString(char * s) const;
  std::string toString() const;
	
  void copyDataInto(void * pointer) const
  {
	// void * memcpy ( void * destination, const void * source, size_t num );
	  memcpy(pointer, longs.data(), sizeof(uint64_t) * (N_LONGS));
  }

  // ABAB: return the raw data packed in an std::array
  // this preserves the lexicographical order on k-mers
  // i.e. A.toString() < B.toString <=> A.getArray() < B.getArray()
  const MERARR &getArray() const {
        return longs;
  }
  const uint8_t *getBytes() const {
        return bytes.data();
  }
  int getNumBytes() const {
        return N_BYTES;
  }
  

  // static functions
  static void set_k(unsigned int _k);
  static constexpr size_t numBytes() { 
	  return (sizeof(uint64_t) * (N_LONGS));
  }


  static const unsigned int MAX_K = MAX_KMER_SIZE;
  static unsigned int k;

 private:
  static unsigned int k_bytes;
  static unsigned int k_modmask; // int?

  // data fields
  union {
    MERARR longs;
    BYTEARR bytes;
  };
  
  unsigned int length;

  // Unions are very useful for low-level programming tasks that involve writing to the same memory area 
  // but at different portions of the allocated memory space, for instance:
  //		union item {
  //			// The item is 16-bits
  //			short theItem;
  //			// In little-endian lo accesses the low 8-bits -
  //			// hi, the upper 8-bits
  //			struct { char lo; char hi; } portions;
  //		};
  //  item tItem;
  //  tItem.theItem = 0xBEAD;
  //  tItem.portions.lo = 0xEF; // The item now equals 0xBEEF
	

 // void shiftForward(int shift);
 // void shiftBackward(int shift);
};


struct KmerHash {
  size_t operator()(const Kmer &km) const {
    return km.hash();
  }
};

// specialization of std::Hash

namespace std
{
    template<> struct hash<Kmer>
    {
        typedef std::size_t result_type;
        result_type operator()(Kmer const& km) const
        {
            auto myhash = km.hash();
            return myhash;
        }
    };
    template<> struct hash<Kmer::MERARR>
    {
        typedef std::size_t result_type;
        result_type operator()(const Kmer::MERARR & km) const
        {
            return MurmurHash3_x64_64((const void*)km.data(),sizeof(Kmer::MERARR));
        }
    };
};


inline std::ostream& operator<<(std::ostream& out, const Kmer& k){
    return out << k.toString();
};

inline std::string toHOPC(std::string original) {
  std::string hopc = "";

  char last = '\n';
  for(auto c : original) {
    if ( last != c ) {
      last = c;
      hopc += c;
    }
  }

  return hopc;
}

inline std::string kHOPC(std::string seq, int len) {
  std::string hopc = "";

  char last = '\n';
  for (auto c : original) {
    if ( last != c ) {
      last = c;
      hopc += c;
    }
    if (hopc.length() == len) return hopc;
  }

  return "";
}

#endif
