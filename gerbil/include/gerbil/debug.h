/*********************************************************************************
Copyright (c) 2016 Marius Erbert, Steffen Rechner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*********************************************************************************/

#ifndef DEBUG_H_
#define DEBUG_H_

#include <iostream>
#include <queue>
#include "config.h"

namespace gerbil {

//#define DEB_CHECK
//#define DEB_MESS

//#define QUEUE_STAT

#ifdef QUEUE_STAT
#define IF_QUEUE_STAT(x) x
#else
#define IF_QUEUE_STAT(x)
#endif

//CLOCK_PROCESS_CPUTIME_ID
//CLOCK_REALTIME
//CLOCK_THREAD_CPUTIME_ID
class StopWatch {
	struct timespec _start;
	struct timespec _end;
	clockid_t _mode;

	__syscall_slong_t _nsec;
	__time_t _sec;
public:
	StopWatch(const clockid_t mode) :
			_mode(mode), _nsec(0), _sec(0) {
	}
	StopWatch() :
			_mode(CLOCK_REALTIME), _nsec(0), _sec(0) {
	}
	void setMode(const clockid_t mode) {
		_mode = mode;
	}

	void start() {
		clock_gettime(_mode, &_start);
		_nsec = 0;
		_sec = 0;
	}
	void stop() {
		clock_gettime(_mode, &_end);
		_nsec += _end.tv_nsec - _start.tv_nsec;
		_sec += _end.tv_sec - _start.tv_sec;
	}
	void hold() {
		clock_gettime(_mode, &_end);
		_nsec += _end.tv_nsec - _start.tv_nsec;
		_sec += _end.tv_sec - _start.tv_sec;
	}
	void proceed() {
		clock_gettime(_mode, &_start);
	}
	double get_s() {
		return _sec + ((double) _nsec / 1000000000);
	}
	double get_ms() {
		return get_s() * 1000;
	}
	double get_us() {
		return get_ms() * 1000;
	}
};

class StackStopWatch {
	std::queue<struct timespec> _starts;
	struct timespec _start;
	struct timespec _end;
	clockid_t _mode;

	bool _inProgress;

	__syscall_slong_t _nsec;
	__time_t _sec;

public:
	StackStopWatch(const clockid_t mode) :
			_mode(mode), _nsec(0), _sec(0), _inProgress(false) {
	}
	StackStopWatch() :
			_mode(CLOCK_REALTIME), _nsec(0), _sec(0), _inProgress(false) {
	}

	void start() {
		clock_gettime(_mode, &_start);
		_starts.push(_start);
		_nsec = 0;
		_sec = 0;
	}
	void stop() {
		clock_gettime(_mode, &_end);
		_start = _starts.front();
		_starts.pop();
		_nsec += _end.tv_nsec - _start.tv_nsec;
		_sec += _end.tv_sec - _start.tv_sec;
	}
	void hold() {
		clock_gettime(_mode, &_end);
		_start = _starts.front();
		_starts.pop();
		_nsec += _end.tv_nsec - _start.tv_nsec;
		_sec += _end.tv_sec - _start.tv_sec;
	}
	void proceed() {
		clock_gettime(_mode, &_start);
		_starts.push(_start);
	}
	double get_s() {
		return _sec + ((double) _nsec / 1000000000);
	}
	double get_ms() {
		return get_s() * 1000;
	}
	double get_us() {
		return get_ms() * 1000;
	}
};

#ifdef DEB
#define DEB_LOG
#define DEB_CHECK

#define DEB_ONLY(x) x
#else
#define DEB_ONLY(x)
#endif

#ifdef DEB_LOG
#define DEB_SS_LOG
#define DEB_KH_LOG
#endif

#ifdef DEB_CHECK
#define CHECK(cond, message) if(!(cond)) printf("CHECK: %s\n\t\t%s (line %d)\n", message, __FILE__, __LINE__);
#else
#define CHECK(cond, message)
#endif

#ifdef DEB_SS_LOG
#define DEB_SS_LOG_RB
#define DEB_SS_LOG_VAL
#define DEB_SS_LOG_SPLIT
#define DEB_SS_LOG_SMER
#endif

#ifdef DEB_MESS
#define DEB_MESS_FASTREADER
#define DEB_MESS_FASTPARSER
#define DEB_MESS_SEQUENCESPLITTER
#define DEB_MESS_SUPERWRITER
#define DEB_MESS_SUPERREADER
#define DEB_MESS_SUPERSPLITTER
#define DEB_MESS_HASHER
#define DEB_MESS_KMCWRITER
#endif

#ifdef DEB_MESS_FASTPARSER
#define IF_MESS_FASTPARSER(x) x
#else
#define IF_MESS_FASTPARSER(x)
#endif

#ifdef DEB_MESS_FASTREADER
#define IF_MESS_FASTREADER(x) x
#else
#define IF_MESS_FASTREADER(x)
#endif

#ifdef DEB_MESS_SEQUENCESPLITTER
#define IF_MESS_SEQUENCESPLITTER(x) x
#else
#define IF_MESS_SEQUENCESPLITTER(x)
#endif

#ifdef DEB_MESS_SUPERWRITER
#define IF_MESS_SUPERWRITER(x) x
#else
#define IF_MESS_SUPERWRITER(x)
#endif

#ifdef DEB_MESS_SUPERREADER
#define IF_MESS_SUPERREADER(x) x
#else
#define IF_MESS_SUPERREADER(x)
#endif

#ifdef DEB_MESS_SUPERSPLITTER
#define IF_MESS_SUPERSPLITTER(x) x
#else
#define IF_MESS_SUPERSPLITTER(x)
#endif

#ifdef DEB_MESS_HASHER
#define IF_MESS_HASHER(x) x
#else
#define IF_MESS_HASHER(x)
#endif

#ifdef DEB_MESS_KMCWRITER
#define IF_MESS_KMCWRITER(x) x
#else
#define IF_MESS_KMCWRITER(x)
#endif

void printChars(char* const a, const uint32_t &l);

void printCharsN(char* const a, const uint32_t &l);

void printByteCodedSeq(const unsigned char* a, const unsigned int &l);

void printByteCodedSeqN(unsigned char* a, const unsigned int &l);

void printByteCodedSeqNT(unsigned char* a, const unsigned int &l,
		const unsigned int &t);

char* getByteCodedSeq(const unsigned char* a, const unsigned int &l);

char* getInt32CodedSeq(const unsigned int &a, const unsigned int &l);
void printInt32CodedSeq(const unsigned int &a, const unsigned int &l);

#define DEB_startThreadClock() struct timespec start_thread, end_thread;\
clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start_thread) ;

#define DEB_stopThreadClock(s) clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_thread) ; \
double dtime_thread = (end_thread.tv_sec - start_thread.tv_sec) + ((double)(end_thread.tv_nsec - start_thread.tv_nsec) / 1000000000); \
printf("Time_%s: %16.6fs\n", s, dtime_thread);

#define DEB_stopThreadIDClock(s, id) clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end_thread) ; \
double dtime_thread = (end_thread.tv_sec - start_thread.tv_sec) + ((double)(end_thread.tv_nsec - start_thread.tv_nsec) / 1000000000); \
printf("Time_%s[%2d]: %16.6fs\n", s, id, dtime_thread);

}

#endif /* DEBUG_H_ */
