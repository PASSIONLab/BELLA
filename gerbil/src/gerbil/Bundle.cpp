/*
 * Bundle.cpp
 *
 *  Created on: 20.05.2015
 *      Author: marius
 */

#include "../../include/gerbil/Bundle.h"

uint gerbil::ReadBundle::K = 0;

gerbil::ReadBundle::ReadBundle(){
	data = new byte[READ_BUNDLE_SIZE_B];
	readsCount = (uint32*)(data + READ_BUNDLE_SIZE_B - 4);
	readOffsets = (uint32*)(data + READ_BUNDLE_SIZE_B - 8);
	*readsCount = 0;
	readOffsets[0] = 0;
}

gerbil::ReadBundle::~ReadBundle(){
	delete[] data;
}

bool gerbil::ReadBundle::isEmpty() const{
	return !*readsCount;
}

void gerbil::ReadBundle::clear() {
	*readsCount = 0;
}

bool gerbil::ReadBundle::add(const uint32 &length, char* read){
	uint32 rc= *readsCount;
	uint32 offset = *(readOffsets - rc);
	if(offset + length + 1 >= READ_BUNDLE_SIZE_B - 8 - rc * 4 - 4)
		return false;
#ifdef USE_MEM_COPY
	memcpy(read, data + offset, length);
#else
	std::copy(read, read + length, data + offset);
#endif
	data[offset + length] = 'E';
	*readsCount = ++rc;
	*(readOffsets-rc) = offset + length + 1;
	return true;
}


bool gerbil::ReadBundle::expand(const uint32 &expLength, char* expRead) {
	uint32 rc= *readsCount;
	uint32 offset = *(readOffsets-rc);
	if(offset + expLength >= READ_BUNDLE_SIZE_B - 8 - rc * 4 - 4)
		return false;
#ifdef USE_MEM_COPY
	memcpy(expRead, data + offset - 1, expLength);
#else
	std::copy(expRead, expRead + expLength, data + offset - 1);
#endif
	*(readOffsets-rc) += expLength;
	data[*(readOffsets-rc) - 1] = 'E';
	return true;
}

bool gerbil::ReadBundle::transfer(ReadBundle* readbundle) {
	uint32 rc= *readsCount;
	bool res = readbundle->add(*(readOffsets - rc) - *(readOffsets - (rc - 1)) - 1, (char*)data + *(readOffsets - (rc - 1)));
	if(!res)
		return false;
	--(*readsCount);
	return true;
}

bool gerbil::ReadBundle::transferKm1(ReadBundle *readbundle) {
	uint32 rc = *readsCount;
	if (*(readOffsets - rc) - *(readOffsets - (rc - 1)) < K) {
		return transfer(readbundle);
	}
	return readbundle->add(K - 1, (char*)data + *(readOffsets - rc) - K);
}

//IF_DEB_DEV(
	void gerbil::ReadBundle::print() {
		uint32 rc= *readsCount;
		printf("readsCount: %6d\n", rc);
		printf("filled: %6u B of %6lu B\n", *(readOffsets - (*readsCount)) + (2 + *readsCount) * 4, READ_BUNDLE_SIZE_B);
		for(uint32 i = 0; i < rc; i++) {
			uint32 l = *(readOffsets - (i + 1)) - *(readOffsets - i);
			printf("[%3d]:  length = %6d (+1)\n", i, l-1);
			char* a = new char[l + 1];
	#ifdef USE_MEM_COPY
			memcpy(data + *(readOffsets - i), a, *(readOffsets - (i + 1)));
	#else
			std::copy(data + *(readOffsets - i), data + *(readOffsets - (i + 1)), a);
	#endif
			a[l] = '\0';
			std::cout << "\t" << a << std::endl;
			delete[] a;
		}

	}
//)

gerbil::SuperBundle::SuperBundle() {
	clear();
}

gerbil::SuperBundle::~SuperBundle() {
}

void gerbil::SuperBundle::finalize() {
	*_next = 0;
	_finalized = true;
}

bool gerbil::SuperBundle::merge(const SuperBundle &sb) {
	if(_next - data + sb._next - sb.data >= SUPER_BUNDLE_DATA_SIZE_B)
		return false;
	_next = std::copy((byte*)sb.data, sb._next, _next);
	sMerNumber += sb.sMerNumber;
	kMerNumber += sb.kMerNumber;
	return true;
}

void gerbil::SuperBundle::clear() {
	_next = data;
	*data = 0;
	_finalized = false;
	tempFileId = TEMPFILEID_NONE;
	tempFileRun = 0;

	sMerNumber = 0;
	kMerNumber = 0;
}


//super-mer in bytes (4x base ==> 1 byte)
bool gerbil::SuperBundle::add(const byte* smer, const uint16 &length, const uint32_t &k) {
	//CHECK(0 < length && length < 128 * 256, "invalid length");
	if(_next - data + (length >> 2) >= SUPER_BUNDLE_DATA_SIZE_B - 4)
		return false;
	++sMerNumber;
	kMerNumber += length - k + 1;
	if(length < 128)
			*(_next++) = length;
	else {
		*(_next++) = (length >> 8) | 0x80;
		*(_next++) = length & 0xff;
	}
	const gerbil::byte* smer_p_end = smer + (length & 0xfffc);
	byte* smer_p = (byte*) smer;

	while(smer_p < smer_p_end) {
		*_next = *smer_p << 6;
		*_next |= *(++smer_p) << 4;
		*_next |= *(++smer_p) << 2;
		*(_next++) |= *((++smer_p)++);
	}
	switch(length & 0x3) {
	case 3: *(_next++) = (*smer_p << 6) | (*(++smer_p) << 4) | (*(++smer_p) << 2); break;
	case 2: *(_next++) = (*smer_p << 6) | (*(++smer_p) << 4); break;
	case 1: *(_next++) = (*smer_p << 6); break;
	}
	return true;
}


bool gerbil::SuperBundle::next(byte* &smer, uint16 &length) {
	if(!(length = *_next))
		return false;
	if(length >= 128)
		length = ((length & 0x7f) << 8) + *(++_next);
	smer = ++_next;
	_next += ((length-1) >> 2) + 1;
	return true;
}


bool gerbil::SuperBundle::isEmpty() const{
	return !*data;
}




gerbil::KmcBundle::KmcBundle() {
	_data = new byte[KMC_BUNDLE_DATA_SIZE_B];
	clear();
}

gerbil::KmcBundle::~KmcBundle() {
	delete[] _data;
}

void gerbil::KmcBundle::clear() {
	_next = _data;
}
