//
// Created by rechner on 5/9/17.
//

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

using namespace std;

// type definitions
typedef unsigned char byte;
typedef unsigned int uint32;

// prints some info
void printHelp() {
	cout << "\ttoFasta <input> <k> [<output>]" << endl;
}

int main(int argc, char** argv) {
	// check the number or parameters
	if(argc < 3 || argc > 4) {
		printHelp();
		return 1;
	}

	// get the parameters
	const string inFilePath(argv[1]);
	const uint32 k = stoi(argv[2]);
	const string outFilePath(argc == 4 ? argv[3] : "");

	// number of bytes per kmer
	const size_t kMerSize_B = (k + 3) / 4;

	// convert two bits to a single base
	const char c[4] = {'A', 'C', 'G', 'T'};

	// the k-mer string
	char kmerSeq[k + 1]; kmerSeq[k] = '\0';

	// files
	FILE* inFile  = NULL;
	FILE* outFile = NULL;

	// open input file
	inFile = fopen ( inFilePath.c_str(), "rb" );

	// check the successful opening of the file
	if(!inFile){
		cerr << "ERROR: Opening the input file '" << inFilePath << "' failed" << endl;
		return 2;
	}

	// get size of input file
	fseek(inFile, 0, SEEK_END);
	const size_t fileSize = ftell(inFile);
	fseek(inFile, 0, SEEK_SET);

	// remaining bytes, which have not been read
	size_t remainingBytes = fileSize;

	if(argc == 4) {
		// open output file
		outFile = fopen ( outFilePath.c_str(), "w+" );

		// check the successful creation of the file
		if(!outFile){
			cerr << "Error: The creation of the output file '" << outFilePath << "' has failed" << endl;
			return 23;
		}
	}

	// print a summary
	printf("input file  : '%s' (%lu B)\n", inFilePath.c_str(), fileSize);
	if(outFile)
		printf("output file : '%s'\n", outFilePath.c_str());
	printf("k           : %u\n", k);
	printf("start converting to FASTA...\n");

	// minimum number of remaining bytes of the buffer (additionally for overlaps)
	size_t minRemBufferBytes = kMerSize_B + 4;

	const size_t bufferSize = 32 * 1024;          // size of the buffer
	size_t bufferOffset     = 0;                  // buffer offset (remaining unused bytes from last turn)
	byte buffer[bufferSize + minRemBufferBytes];  // buffer
	size_t curSize;                               // number of currently read bytes
	size_t i                = 0;                  // progress purpose

	// disable console buffering
	setbuf(stdout, NULL);

	// main loop
	while(remainingBytes) {
		// progress
		if(outFile && !(i++ % 1024))
			printf("\r%lu B left", remainingBytes);

		// determine number of readable bytes
		if(remainingBytes > bufferSize) {
			curSize = bufferSize;
			remainingBytes -= bufferSize;
		}
		else {
			curSize = remainingBytes;
			remainingBytes = 0;
		}
		// read bytes from file into the buffer
		if(fread(buffer + bufferOffset, 1, curSize, inFile) != curSize) {
			cerr << "ERROR: An error has occurred while reading bytes from the input file!" << endl;
			return 3;
		}

		// initialize iterator and final pointer
		const byte* p   = buffer;
		const byte* end = buffer + curSize + bufferOffset;

		// k-mer counter
		uint32 counter;

		// loop over buffer
		while(remainingBytes ? p + minRemBufferBytes <= end : p < end) {
			// get counter value (small)
			counter = (uint32)*(p++);
			if(counter >= 255) {
				// large value
				counter = *((uint32*)p);
				p += 4;
			}
			// k-mer: convert bytes to string
			for(uint i = 0; i < k; ++i)
				kmerSeq[i] = c[(p[i >> 2] >> (2 * (3 - (i & 0x3)))) & 0x3];

			// increase pointer
			p += kMerSize_B;

			// print fasta (console/file)
			if(outFile)
				fprintf(outFile, ">%u\n%s\n", counter, kmerSeq);
			else
				printf(">%u\n%s\n", counter, kmerSeq);
		}

		// carryover (remaining bytes: possibly incomplete k-mers with their counters)
		if(p < end) {
			bufferOffset = end - p;
			// copy bytes to the beginning
			std::memcpy(buffer, p, bufferOffset);
		}
		else
			bufferOffset = 0;
	}

	// print read/written bytes
	printf("\rbytes read    : %lu B          \n", fileSize);
	printf("\rbytes written : %lu B          \n", ftell(outFile));

	// close files
	fclose(inFile);
	if(outFile) fclose(outFile);

	// exit without errors
	return 0;
}