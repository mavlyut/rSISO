#ifndef DEFINES
#define DEFINES

#include <chrono>
#include <fstream>
#include <random>

#define ENABLE_OPT_1
#define ENABLE_OPT_2
#define ENABLE_OPT_3
#define ENABLE_OPT_5

// #define LOG
// #define TIMELOG
// #define CNTLOG
// #define TEST
#define EPS 1e-10

static const _Float64 INF = INFINITY;
static constexpr _Float64 COEF = M_2_SQRTPI * M_SQRT1_2 / 2;
inline std::mt19937 gen(time(0));

static const unsigned UNINIT = -1;
static const unsigned chunk_size = 64;

typedef std::chrono::milliseconds ms;

#if defined(LOG) && !defined(TEST)
inline std::ofstream _log_out("_log_out");
#define __log(msg) _log_out << msg;
#else
#define __log(msg)
#endif

#define fail(b, msg)							\
	if (!(b)) { 								\
		__log("FAILURE: " << msg << std::endl);	\
	}

static std::size_t SUM_CNT;
static std::size_t MUL_CNT;
static std::size_t CNT_BIN;

void clear_cnt();
void clear_cnt_bin();

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<T> const& a) {
	for (unsigned i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << ' ';
		}
		out << a[i];
	}
	return out;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, std::vector<std::vector<T>> const& a) {
	for (unsigned i = 0; i < a.size(); i++) {
		if (i != 0) {
			out << "\n\t";
		}
		out << a[i];
	}
	return out;
}

template <typename T>
std::istream& operator>>(std::istream& in, std::vector<T>& a) {
	for (T& i : a) {
		in >> i;
	}
	return in;
}

#endif // DEFINES
