#ifndef GLOBAL_H_
#define GLOBAL_H_

struct readType_ {

	std::string nametag;   
	std::string seq; 
	int readid;

    bool operator < (readType_ & str)
    {
        return (readid < str.readid);
    }

};

typedef vector<readType_> readVector_;

#ifdef __cplusplus
extern "C" {
#endif

#include "../optlist/optlist.h" // command line parser
#include "../libgaba/gaba.h" // sequence alignment lib

#ifdef __cplusplus
}
#endif

#include "../libcuckoo/cuckoohash_map.hh"

#endif