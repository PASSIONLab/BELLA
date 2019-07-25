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

#ifndef SYNCQUEUE_H_
#define SYNCQUEUE_H_

#include "types.h"
#include <queue>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <condition_variable>
#include <boost/thread/thread.hpp>
#include <atomic>
#include <cstddef>

namespace gerbil {

// base class for special synchronized queues
class SyncSwapQueue {
protected:
	bool _isFinalized;

	std::atomic<size_t>  _tail;
	std::atomic<size_t>  _head;

	const uint64 _size;
	uint64 _mask;

	SyncSwapQueue(const uint64 size):_isFinalized(false), _tail(0), _head(0), _size(size), _mask(size - 1){
		assert(_size > 1);
		// 2^x size only
		//assert((_size & _mask) == 0);
		assert(_tail.is_lock_free());
		assert(_head.is_lock_free());
	};
	~SyncSwapQueue() {
	};

	inline size_t increment(const size_t& idx) const {
	  return (idx + 1) % _size; // & _mask
	}

public:
	inline bool empty() const{
	  return _head.load() == _tail.load();
	}

	//prod only
	inline void finalize() {
		_isFinalized = true;
	}

	inline void reset() {
		assert(_isFinalized);
		_isFinalized = false;
	}

	inline const bool& isFinalized() const {
		return _isFinalized;
	}
};


#define SSQ_SPSC_WAIT sched_yield()
//#define SSQ_SPSC_WAIT usleep(1);
template <typename T>
class SyncSwapQueueSPSC: public SyncSwapQueue{
	T** _data;
	IF_QUEUE_STAT(
		uint64 _wPush;
		uint64 _wPop;
		uint64 _ops;
		StopWatch _swPop;
		StopWatch _swPush;
		uint64 _maxUsedSize;
	)
public:
	SyncSwapQueueSPSC(const uint64 size)
		:SyncSwapQueue(size)
	{
		_data = new T*[_size];
		T** p_end(_data + _size);
		for(T** p(_data); p < p_end; ++p)
			*p = new T;
		IF_QUEUE_STAT(
			_wPush = _wPop = _ops = 0;
			_maxUsedSize = 0;
			_swPop.start();
			_swPop.hold();
			_swPush.start();
			_swPush.hold();
		)
	}

	~SyncSwapQueueSPSC(){
		T** p_end = _data + _size;
		for(T** p(_data); p < p_end; ++p)
			delete *p;
		delete[] _data;
		IF_QUEUE_STAT(
			_swPop.proceed();
			_swPush.proceed();
			_swPop.stop();
			_swPush.stop();
			printf("SPSC:   ");
			printf("max: %5lu/%5lu\t", _maxUsedSize, _size);
			printf("push: %.3f s   pop: %.3f s    | ", _swPush.get_s(), _swPop.get_s());
			printf("wPush=%9lu  wPop=%9lu  ops=%9lu\n", _wPush, _wPop, _ops);
		)
	}

	void swapPush(T* &item) {
		IF_QUEUE_STAT(_swPush.proceed();)
		assert(!_isFinalized);
		const auto current_head = _head.load(std::memory_order_relaxed);
		const size_t next_head = increment(current_head);
		IF_QUEUE_STAT(++_ops;)
		while(next_head == _tail.load(std::memory_order_acquire)){
			IF_QUEUE_STAT(++_wPush;)
			SSQ_SPSC_WAIT;
		}
		std::swap(_data[current_head], item);
		_head.store(next_head, std::memory_order_release);
		IF_QUEUE_STAT(
			uint64 size((_size + _head.load() - _tail.load()) % _size);
			if(_maxUsedSize < size)
				_maxUsedSize = size;
		)
		IF_QUEUE_STAT(_swPush.hold();)
	}


	//tauscht die Elemente aus und rueckt eins weiter
	bool swapPop(T* &item) {
		IF_QUEUE_STAT(_swPop.proceed();)
		const auto current_tail = _tail.load(std::memory_order_relaxed);
		std::atomic<size_t> current_head;
		while(current_tail == (current_head = _head.load(std::memory_order_acquire)) && !_isFinalized) {
			IF_QUEUE_STAT(++_wPop;)
			//SSQ_SPSC_WAIT;
			usleep(FAST_BLOCK_SIZE_B / 1024 * (MIN_FASTBUNDLEBUFFER_SIZE_B / FAST_BLOCK_SIZE_B / 32));
		}
		if(current_tail == current_head) {
			assert(_isFinalized);
			IF_QUEUE_STAT(_swPop.hold();)
			return false;
		}
		std::swap(_data[current_tail], item);

		_tail.store(increment(current_tail), std::memory_order_release);
		IF_QUEUE_STAT(_swPop.hold();)
		return true;
	}
};


template <typename T>
class SyncSwapQueueMPSC{
	typedef enum {sq_empty, sq_fill} TSQState;

	T** _data;
	T** _prod_d_p;
	T** _cons_d_p;

	byte* _state;
	byte* _prod_s_p;
	byte* _cons_s_p;
	uint64 _maxSize;
	uint64 _size;

	std::mutex _mtx;
	std::condition_variable _cv_empty;
	std::condition_variable _cv_full;
	bool _finalize;

	IF_QUEUE_STAT(
		uint64 _wPush;
		uint64 _wPop;
		uint64 _ops;
		StackStopWatch _swPush;
		StopWatch _swPop;
		uint64 _maxUsedSize;
	)
public:
	SyncSwapQueueMPSC(const uint64 maxSize){
		_finalize = false;
		_maxSize = maxSize;
		_size = 0;
		_data = new T*[_maxSize];
		_prod_d_p = _cons_d_p = _data;
		_state = new byte[_maxSize];
		_prod_s_p = _cons_s_p = _state;
		for(uint i = 0; i < _maxSize; i++) {
			_state[i] = sq_empty;
			_data[i] = new T();
		}

		IF_QUEUE_STAT(
			_wPush = _ops = 0;
			_wPop = 0;
			_maxUsedSize = 0;
			_swPop.start();
			_swPop.hold();
			_swPush.start();
			_swPush.hold();
		)
	}

	~SyncSwapQueueMPSC(){
		for(uint i = 0; i < _maxSize; i++)
			delete _data[i];

		delete[] _state;
		delete[] _data;
		IF_QUEUE_STAT(
			_swPop.proceed();
			_swPush.proceed();
			_swPop.stop();
			_swPush.stop();
			printf("MPSC:   ");
			printf("max: %5lu/%5lu\t", _maxUsedSize, _maxSize);
			printf("push: %.3f s   pop: %.3f s    | ", _swPush.get_s(), _swPop.get_s());
			printf("wPush=%9lu  wPop=%9lu  ops=%9lu\n", _wPush, _wPop, _ops);
		)
	}
	void clear() {
		_finalize = false;
		_prod_d_p = _cons_d_p = _data;
		_prod_s_p = _cons_s_p = _state;
		for(uint i = 0; i < _maxSize; i++)
			_state[i] = sq_empty;
	}

	void reset() {
		_finalize = false;
	}

	inline bool empty() {
	  return !_size;
	}

	void finalize() {
		_finalize = true;
		_cv_empty.notify_one();
	}

	void swapPush(T* &item) {
		if(_finalize)
			return;

		std::unique_lock<std::mutex> ulock(_mtx);
		IF_QUEUE_STAT(
			++_ops;
		)
		while(*_prod_s_p == sq_fill){
			IF_QUEUE_STAT(
				++_wPush;
				_swPush.proceed();
			)
			_cv_full.wait(ulock);
			IF_QUEUE_STAT(
				_swPush.hold();
			)
		}

		T* swap = *_prod_d_p;
		*(_prod_d_p++) = item;
		item = swap;

		*(_prod_s_p++) = sq_fill;
		if(_prod_d_p >= _data + _maxSize) {
			_prod_d_p -= _maxSize;
			_prod_s_p -= _maxSize;
		}
		++_size;
		IF_QUEUE_STAT(
			if(_maxUsedSize < _size)
				_maxUsedSize = _size;
		)
		ulock.unlock();
		_cv_empty.notify_one();
	}

	/**
	 * It is completely unclear what the return value of this method is.
	 */
	bool swapPop(T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		while(!_finalize && !_size) {
			IF_QUEUE_STAT(
				++_wPop;
				_swPop.proceed();
			)
			_cv_empty.wait(ulock);
			IF_QUEUE_STAT(
				_swPop.hold();
			)
		}
		if(!_size)
			return false;
		T* swap = *_cons_d_p;
		*(_cons_d_p++) = item;
		item = swap;
		*(_cons_s_p++) = sq_empty;
		if(_cons_d_p >= _data + _maxSize) {
			_cons_d_p -= _maxSize;
			_cons_s_p -= _maxSize;
		}
		--_size;
		ulock.unlock();
		_cv_full.notify_one();

		return true;
	}

	bool swapPop_nl(T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		if(!_size)
			return false;

		T* swap = *_cons_d_p;
		*(_cons_d_p++) = item;
		item = swap;

		*(_cons_s_p++) = sq_empty;
		if(_cons_d_p >= _data + _maxSize) {
			_cons_d_p -= _maxSize;
			_cons_s_p -= _maxSize;
		}
		--_size;
		ulock.unlock();
		_cv_full.notify_one();

		return true;
	}

	inline bool isFinalized() {
		return _finalize;
	}
};

template <typename T>
class SyncSwapQueueSPMC{
	typedef enum {sq_empty, sq_fill} TSQState;

	T** _data;
	T** _prod_d_p;
	T** _cons_d_p;

	byte* _state;
	byte* _prod_s_p;
	byte* _cons_s_p;
	uint64 _maxSize;
	uint64 _size;

	std::mutex _mtx;
	std::condition_variable _cv_empty;
	std::condition_variable _cv_full;

	bool _finalize;

	IF_QUEUE_STAT(
		uint64 _wPush;
		uint64 _wPop;
		uint64  _ops;
		StopWatch _swPush;
		StackStopWatch _swPop;
		uint64 _maxUsedSize;
	)
public:
	SyncSwapQueueSPMC(const uint64 maxSize){
		_finalize = false;
		_maxSize = maxSize;
		_size = 0;
		_data = new T*[_maxSize];
		_prod_d_p = _cons_d_p = _data;
		_state = new byte[_maxSize];
		_prod_s_p = _cons_s_p = _state;
		for(uint i = 0; i < _maxSize; ++i) {
			_state[i] = sq_empty;
			_data[i] = new T();
		}
		IF_QUEUE_STAT(
			_wPop = _ops = 0;
			_wPush = 0;
			_maxUsedSize = 0;
			_swPop.start();
			_swPop.hold();
			_swPush.start();
			_swPush.hold();
		)
	}

	~SyncSwapQueueSPMC(){
		/*for(_cons_d_p = _data, _prod_d_p = _data + _maxSize; _cons_d_p < _prod_d_p; _cons_d_p++)
			if(*_cons_d_p)
				delete *_cons_d_p;*/
		for(uint i = 0; i < _maxSize; ++i)
			delete _data[i];
		delete[] _state;
		delete[] _data;
		IF_QUEUE_STAT(
			_swPop.proceed();
			_swPush.proceed();
			_swPop.stop();
			_swPush.stop();
			printf("SPMC:   ");
			printf("max: %5lu/%5lu\t", _maxUsedSize, _maxSize);
			printf("push: %.3f s   pop: %.3f s    | ", _swPush.get_s(), _swPop.get_s());
			printf("wPush=%9lu  wPop=%9lu  ops=%9lu\n", _wPush, _wPop, _ops);
		)
	}

	inline bool empty() {
	  return !_size;
	}

	void finalize() {
		_finalize = true;
		_cv_empty.notify_all();
	}

	void swapPush(T* &item) {
		if(_finalize)
			return;

		std::unique_lock<std::mutex> ulock(_mtx);
		IF_QUEUE_STAT(
			++_ops;
		)
		while(*_prod_s_p == sq_fill) {
			IF_QUEUE_STAT(
				++_wPush;
				_swPush.proceed();
			)
			_cv_full.wait(ulock);
			IF_QUEUE_STAT(
				_swPush.hold();
			)
		}

		T* swap = *_prod_d_p;
		*(_prod_d_p++) = item;
		item = swap;

		*(_prod_s_p++) = sq_fill;
		if(_prod_d_p >= _data + _maxSize) {
			_prod_d_p -= _maxSize;
			_prod_s_p -= _maxSize;
		}
		++_size;
		IF_QUEUE_STAT(
			if(_maxUsedSize < _size)
				_maxUsedSize = _size;
		)
		ulock.unlock();
		_cv_empty.notify_one();
	}

	bool swapPop(T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		while(!_finalize && !_size) {
			IF_QUEUE_STAT(
				++_wPop;
				_swPop.proceed();
			)
			_cv_empty.wait(ulock);
			IF_QUEUE_STAT(
				_swPop.hold();
			)
		}
		if(!_size)
			return false;

		T* swap = *_cons_d_p;
		*(_cons_d_p++) = item;
		item = swap;

		*(_cons_s_p++) = sq_empty;
		if(_cons_d_p >= _data + _maxSize) {
			_cons_d_p -= _maxSize;
			_cons_s_p -= _maxSize;
		}
		--_size;
		ulock.unlock();
		_cv_full.notify_one();

		return true;
	}

	//Zugriff auf das i. item, NICHT gleichzeitig mit Konsumenten verwenden
	bool top(const uint64 i, T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		while(!_finalize && _size <= i)
			_cv_empty.wait(ulock);
		if(_size <= i)
			return false;
		T** p = _cons_d_p + i;
		if(p >= _data + _maxSize)
			p -= _maxSize;
		item = *p;
		return true;
	}

	inline bool isFinalized() {
		return _finalize;
	}
};

template <typename T>
class SyncSwapQueueMPMC{
	typedef enum {sq_empty, sq_fill} TSQState;

	T** _data;
	T** _prod_d_p;
	T** _cons_d_p;

	byte* _state;
	byte* _prod_s_p;
	byte* _cons_s_p;
	uint64 _maxSize;
	uint64 _size;

	std::mutex _mtx;
	std::condition_variable _cv_empty;
	std::condition_variable _cv_full;

	bool _finalize;

	IF_QUEUE_STAT(
		uint64 _wPush;
		uint64 _wPop;
		uint64  _ops;
		StackStopWatch _swPush;
		StackStopWatch _swPop;
		uint64 _maxUsedSize;
	)
public:
	SyncSwapQueueMPMC(const uint64 maxSize){
		_finalize = false;
		_maxSize = maxSize;
		_size = 0;
		_data = new T*[_maxSize];
		_prod_d_p = _cons_d_p = _data;
		_state = new byte[_maxSize];
		_prod_s_p = _cons_s_p = _state;
		for(uint i = 0; i < _maxSize; ++i) {
			_state[i] = sq_empty;
			_data[i] = new T();
		}
		IF_QUEUE_STAT(
			_wPop = _ops = 0;
			_wPush = 0;
			_maxUsedSize = 0;
			_swPop.start();
			_swPop.hold();
			_swPush.start();
			_swPush.hold();
		)
	}

	~SyncSwapQueueMPMC(){
		/*for(_cons_d_p = _data, _prod_d_p = _data + _maxSize; _cons_d_p < _prod_d_p; _cons_d_p++)
			if(*_cons_d_p)
				delete *_cons_d_p;*/
		for(uint i = 0; i < _maxSize; ++i)
			delete _data[i];
		delete[] _state;
		delete[] _data;
		IF_QUEUE_STAT(
			_swPop.proceed();
			_swPush.proceed();
			_swPop.stop();
			_swPush.stop();
			printf("MPMC:   ");
			printf("max: %5lu/%5lu\t", _maxUsedSize, _maxSize);
			printf("push: %.3f s   pop: %.3f s    | ", _swPush.get_s(), _swPop.get_s());
			printf("wPush=%9lu  wPop=%9lu  ops=%9lu\n", _wPush, _wPop, _ops);
		)
	}

	inline bool empty() {
	  return !_size;
	}

	void finalize() {
		_finalize = true;
		_cv_empty.notify_all();
	}

	void swapPush(T* &item) {
		if(_finalize)
			return;

		std::unique_lock<std::mutex> ulock(_mtx);
		IF_QUEUE_STAT(
			++_ops;
		)
		while(*_prod_s_p == sq_fill) {
			IF_QUEUE_STAT(
				++_wPush;
				_swPush.proceed();
			)
			_cv_full.wait(ulock);
			IF_QUEUE_STAT(
				_swPush.hold();
			)
		}

		T* swap = *_prod_d_p;
		*(_prod_d_p++) = item;
		item = swap;

		*(_prod_s_p++) = sq_fill;
		if(_prod_d_p >= _data + _maxSize) {
			_prod_d_p -= _maxSize;
			_prod_s_p -= _maxSize;
		}
		++_size;
		IF_QUEUE_STAT(
			if(_maxUsedSize < _size)
				_maxUsedSize = _size;
		)
		ulock.unlock();
		_cv_empty.notify_one();
	}

	bool swapPop(T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		while(!_finalize && !_size) {
			IF_QUEUE_STAT(
				++_wPop;
				_swPop.proceed();
			)
			_cv_empty.wait(ulock);
			IF_QUEUE_STAT(
				_swPop.hold();
			)
		}
		if(!_size)
			return false;

		T* swap = *_cons_d_p;
		*(_cons_d_p++) = item;
		item = swap;

		*(_cons_s_p++) = sq_empty;
		if(_cons_d_p >= _data + _maxSize) {
			_cons_d_p -= _maxSize;
			_cons_s_p -= _maxSize;
		}
		--_size;
		ulock.unlock();
		_cv_full.notify_one();

		return true;
	}

	//Zugriff auf das i. item, NICHT gleichzeitig mit Konsumenten verwenden
	bool top(const uint64 i, T* &item) {
		std::unique_lock<std::mutex> ulock(_mtx);
		while(!_finalize && _size <= i)
			_cv_empty.wait(ulock);
		if(_size <= i)
			return false;
		T** p = _cons_d_p + i;
		if(p >= _data + _maxSize)
			p -= _maxSize;
		item = *p;
		return true;
	}

	inline bool isFinalized() {
		return _finalize;
	}

	inline const uint64& getMaxSize() {
		return _maxSize;
	}
};

}

#endif /* SYNCQUEUE_H_ */
