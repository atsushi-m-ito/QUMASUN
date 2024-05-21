#pragma once
#include <chrono>

template <typename DURATION>
inline
double DoubleSec(DURATION d) {
	return 1.0e-9 * (double)std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
}

