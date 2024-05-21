


/*
経過時間計測クラスv2

*/
#pragma once

#include <cstdlib>
#include <cstdio>
#include <chrono>

template <int NUM_SECTOR, bool IS_MEASURE>
class StopWatch
{
	
private:
	
	double m_duration[NUM_SECTOR]{ 0.0 };
	int64_t m_count[NUM_SECTOR]{ 0 };
	
	std::chrono::high_resolution_clock::time_point m_prev_tm;


public:

	StopWatch() {		
		m_prev_tm = std::chrono::high_resolution_clock::now();
	}

	~StopWatch() = default;

	double Record(int index, int64_t count = 1) {
		const auto current_tm = std::chrono::high_resolution_clock::now();
		const double sector_tm = DoubleSec(current_tm - m_prev_tm);
		m_duration[index] += sector_tm;
		m_count[index] += count;
		m_prev_tm = current_tm;
		return sector_tm;
	};

	void Restart() {
		m_prev_tm = std::chrono::high_resolution_clock::now();
	};

	void Print(const char* title, int index) const {
		FPrint(stdout, title, index);
	}

	void FPrint(FILE* fp, const char* title, int index) const {

		fprintf(fp, "%s %f [s] / %zd [call]\n", title, m_duration[index], m_count[index]);
	}

	double GetTime(int index) {
		if ((0 <= index) && (index < NUM_SECTOR)) {
			return m_duration[index];
		}
		return 0.0;
	}

	int64_t GetCount(int index) {
		if ((0 <= index) && (index < NUM_SECTOR)) {
			return m_count[index];
		}
		return 0.0;
	}


	double Total() const {
		double sum = 0.0;
		for (int i = 0; i < NUM_SECTOR; ++i) {
			sum += m_duration[i];
		}
		return sum;
	}
	
	double Total(const std::vector<int>& ids) const {
		double sum = 0.0;
		for (const auto& i : ids ){
			sum += m_duration[i];
		}
		return sum;
	}

	void Clear() {
		for (int i = 0; i < NUM_SECTOR; ++i) {
			m_duration[i] = 0.0;
			m_count[i] = 0;
		}
	}

private:

	template <typename DURATION>
	static 
	double DoubleSec(DURATION d) {
		return 1.0e-9 * (double)std::chrono::duration_cast<std::chrono::nanoseconds>(d).count();
	}


};

template <int NUM_SECTOR>
class StopWatch< NUM_SECTOR, false>
{

public:

	StopWatch() = default; 

	~StopWatch() = default;

	double Record(int index) {
		return 0.0;
	};

	void Restart() {};

	void Print(const char* title, int index) const {}

	void FPrint(FILE* fp, const char* title, int index) const {}

	double Total() const {	return 0.0;	}

	void Clear() {}
};



