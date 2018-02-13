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
#endif