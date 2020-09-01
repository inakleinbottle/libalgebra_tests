#include <benchmark/benchmark.h>

#include "basis.h"
#include <random>



typedef alg::vectors::vector<Basis<2<<13>, true> Vec;


Vec generate_combined_vec(size_t num, size_t max_diff) {

    std::mt19937 key_rng;
    std::mt19937 sca_rng;

    std::uniform_real_distribution<double> reals (-2.0, 2.0);

    Vec nvec {};
    size_t crr = 0;
    size_t rn;
    for (size_t i=0; i<num; ++i) {
        //rn = crr + key_rng()%max_diff;
        //crr = rn;
        nvec.add_scal_prod(typename Vec::KEY(i), reals(sca_rng));
    }

    return nvec;
}



static void time_addition_combined_vectors(benchmark::State &state)
{
    Vec vec1 = generate_combined_vec(state.range(0), 25);
    Vec vec2 = generate_combined_vec(state.range(0), 25);
    for (auto _ : state) {
        benchmark::DoNotOptimize(vec1+vec2);
    }
}


BENCHMARK(time_addition_combined_vectors)->Range(2<<6, 2<<15);