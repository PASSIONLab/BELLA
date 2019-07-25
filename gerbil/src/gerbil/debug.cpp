/*
 * debug.cpp
 *
 *  Created on: 20.05.2015
 *      Author: marius
 */

#include "../../include/gerbil/debug.h"

void gerbil::printChars(char* const a, const uint32_t &l) {
	char s[l+1];
	std::copy(a, a+l, s);
	s[l] = '\0';
	std::printf("%s", s);
}

void gerbil::printCharsN(char* const a, const uint32_t &l) {
	char s[l+1];
	std::copy(a, a+l, s);
	s[l] = '\0';
	std::printf("%s\n", s);
}

void gerbil::printByteCodedSeq(const unsigned char* a, const unsigned int &l) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char s[l+1];
	s[l] = '\0';
	for(uint i = 0; i < l; i++)
		s[i] = c[a[i >> 2] & 0x3];
	std::printf("%s", s);
}

void gerbil::printByteCodedSeqN(unsigned char* a, const unsigned int &l) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char s[l+1];
	s[l] = '\0';
	for(uint i = 0; i < l; i++)
		s[i] = c[(a[i >> 2] >> (2 * (3 - (i & 0x3)))) & 0x3];
	std::printf("%s\n", s);
}

char* gerbil::getByteCodedSeq(const unsigned char* a, const unsigned int &l) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char* s = new char[l+1];
	s[l] = '\0';
	for(uint i = 0; i < l; i++)
		s[i] = c[(a[i >> 2] >> (2 * (3 - (i & 0x3)))) & 0x3];
	return s;
}

char* gerbil::getInt32CodedSeq(const unsigned int &a, const unsigned int &l) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char* s = new char[l+1];
	s[l] = '\0';
	for(uint i(0); i <l; i++)
		s[i] = c[(a >> (2 * (l - i - 1))) & 0x3];
	return s;
}

void gerbil::printInt32CodedSeq(const unsigned int &a, const unsigned int &l) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char s[l+1];
	s[l] = '\0';
	for(uint i(0); i <l; i++)
		s[i] = c[(a >> (2 * (l - i - 1))) & 0x3];
	std::printf("%s", s);
}


void gerbil::printByteCodedSeqNT(unsigned char* a, const unsigned int &l, const unsigned int &t) {
	char c[4] = {'A', 'C', 'G', 'T'};
	char s[l+1+t];
	s[l+t] = '\0';
	for(uint i = 0; i < t; i++)
		s[i] = ' ';
	for(uint i = 0; i < l; i++)
		s[i+t] = c[(a[i >> 2] >> (2 * (3 - (i & 0x3)))) & 0x3];
	std::printf("%s\n", s);
}
