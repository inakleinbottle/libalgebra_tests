#include <benchmark/benchmark.h>

#include <libalgebra/libalgebra.h>
#include <random>

#ifndef USE_FLATMAP
#define MAP_T std::map<KEY, SCA>
#else
#include <boost/container/flat_map.hpp>
#define MAP_T boost::container::flat_map<KEY, SCA>
#endif


class Basis {
public:
    typedef double RATIONAL;
    typedef size_t KEY;
    typedef double SCA;
    typedef MAP_T  MAP;

    // Default constructor
    Basis() {}

    unsigned degree(const KEY &k) {
        return 1;
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<Basis*, KEY> &t
    ) {
        return os << t.second;
    }

};
typedef alg::sparse_vector<Basis> Vec;


Vec generate_vec(size_t num, size_t max_diff) {

    std::mt19937 key_rng;
    std::mt19937 sca_rng;

    std::uniform_real_distribution<double> reals (-2.0, 2.0);

    Vec nvec {};
    size_t crr = 0;
    size_t rn;
    for (size_t i=0; i<num; ++i) {
        rn = crr + key_rng()%max_diff;
        crr = rn;
        nvec[crr] = reals(sca_rng);
    }

    return nvec;
}



static void time_addition_vectors(benchmark::State &state)
{
    Vec vec1 = generate_vec(state.range(0), 25);
    Vec vec2 = generate_vec(state.range(0), 25);
    for (auto _ : state) {
        benchmark::DoNotOptimize(vec1+vec2);
    }
}


BENCHMARK(time_addition_vectors)->Range(2<<6, 2<<15);