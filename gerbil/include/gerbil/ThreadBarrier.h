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

#ifndef THREADBARRIER_H_
#define THREADBARRIER_H_

#include <mutex>
#include <condition_variable>

namespace gerbil {

/**
 * @brief Represents a CPU thread barrier
 * @note The barrier automatically resets after all threads are synced
 */
class Barrier {
private:
	std::mutex m_mutex;
	std::condition_variable m_cv;

	size_t m_count;
	const size_t m_initial;

	enum State
		: unsigned char {
			Up, Down
	};
	State m_state;

public:
	explicit Barrier(std::size_t count) :
			m_count { count }, m_initial { count }, m_state { State::Down } {
	}

	/// Blocks until all N threads reach here
	void sync() {
		std::unique_lock < std::mutex > lock { m_mutex };

		if (m_state == State::Down) {
			// Counting down the number of syncing threads
			if (--m_count == 0) {
				m_state = State::Up;
				m_cv.notify_all();
			} else {
				m_cv.wait(lock, [this] {return m_state == State::Up;});
			}
		}

		else // (m_state == State::Up)
		{
			// Counting back up for Auto reset
			if (++m_count == m_initial) {
				m_state = State::Down;
				m_cv.notify_all();
			} else {
				m_cv.wait(lock, [this] {return m_state == State::Down;});
			}
		}
	}
};

}

#endif /* THREADBARRIER_H_ */
