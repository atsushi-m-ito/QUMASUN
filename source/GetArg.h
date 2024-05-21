#pragma once
#include <cstring>
#include <type_traits>

template <class T>
inline
T GetArgumentNum(int argc, char* const argv[], const char* key, const T& default_val) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], key)==0) {
			if (i + 1 >= argc) break;
			if constexpr (std::is_same_v<T, double>) {
				return strtod(argv[i + 1], nullptr);
			}
			if constexpr (std::is_same_v<T, float>) {
				return (float)strtod(argv[i + 1], nullptr);
			}			
			return strtol(argv[i + 1], nullptr, 10);			
		}
	}
	return default_val;
}

template <class T>
inline
int GetArgumentNumList(int argc, char* const argv[], const char* key, int num, T* results) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], key) == 0) {
			int k;
			for (k = 0; k < std::min(num, argc - i - 1); ++k) {
				if constexpr (std::is_same_v < T, double> ) {
					results[k] = strtod(argv[k + i + 1], nullptr);
				}else if constexpr (std::is_same_v<T, float>) {
					results[k] = (float)strtod(argv[k + i + 1], nullptr);
				}
				else {
					results[k] = strtol(argv[k+i + 1], nullptr, 10);
				}
			}

			return k;
		}
	}
	return 0;
}


inline
bool GetArgumentText(int argc, char* const argv[], const char* key, char* res) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], key) == 0) {
			if (i + 1 >= argc) {
				return true;
			}
			strcpy(res, argv[i + 1]);
			return true;
		}
	}
	return false;
}



