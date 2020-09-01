
#include <benchmark/benchmark.h>

#include <libalgebra/libalgebra.h>
#include <random>

typedef alg::free_tensor_basis<double, double, 2, 32> BASIS;
typedef typename BASIS::KEY KEY;

typedef alg::vectors::dense_vector<BASIS> DVec;
typedef alg::vectors::sparse_vector<BASIS> SVec;


template <typename V>
void generate_full_vector(V& vec, size_t dimension)
{
    BASIS b {};
    KEY key{b.begin()};

    std::uniform_real_distribution<double> real(-1.0, 1.0);

    std::mt19937 rng {};

    for (size_t i=0; i<dimension; ++i) {
        vec[key] = real(rng);
        key = b.nextkey(key);
    }
}


static void time_iteration_sparse_vector(benchmark::State &state)
{
    auto dimension = state.range(0);
    SVec vec {};
    generate_full_vector(vec, dimension);
    bool v{false};

    for (auto _ : state) {
        for (auto& item : vec)
            benchmark::DoNotOptimize(v |= (item.second > 0));
    }
    state.SetComplexityN(dimension);
}

static void time_iteration_dense_vector(benchmark::State& state)
{
    auto dimension = state.range(0);
    DVec vec(dimension);
    generate_full_vector(vec, dimension);
    bool v{false};

    for (auto _ : state) {
        for (auto item : vec.get_coeffs())
            benchmark::DoNotOptimize(v |= (item > 0));
    }
    state.SetComplexityN(dimension);
}

BENCHMARK(time_iteration_sparse_vector)
    ->RangeMultiplier(2)
    ->Range(32, 2<<25)
    ->Complexity(benchmark::oN);
BENCHMARK(time_iteration_dense_vector)
    ->RangeMultiplier(2)
    ->Range(32, 2<<25)
    ->Complexity(benchmark::oN);



template <typename V>
void generate_sparse_vector(V& vec, size_t dimension)
{
    BASIS b {};
    KEY key{b.begin()};

    std::uniform_real_distribution<double> real(-1.0, 1.0);

    std::mt19937 rng {dimension};
    std::mt19937 pos_rng {dimension % 12345};
    long long filter = 5 + pos_rng() % 25;

    for (size_t i=0; i<dimension; ++i) {
        if (i % filter == 0)
            vec[key] = real(rng);
        key = b.nextkey(key);
    }
}

static void time_sparse_iteration_sparse_vector(benchmark::State &state)
{
    auto dimension = state.range(0);
    SVec vec {};
    generate_sparse_vector(vec, dimension);
    bool v{false};

    for (auto _ : state) {
        for (auto& item : vec)
            benchmark::DoNotOptimize(v |= (item.second > 0));
    }
    state.SetComplexityN(dimension);
}

static void time_sparse_iteration_dense_vector(benchmark::State& state)
{
    auto dimension = state.range(0);
    DVec vec(dimension);
    generate_sparse_vector(vec, dimension);
    bool v{false};

    for (auto _ : state) {
        for (auto item : vec.get_coeffs())
            benchmark::DoNotOptimize(v |= (item > 0));
    }
    state.SetComplexityN(dimension);
}

BENCHMARK(time_sparse_iteration_sparse_vector)
    ->RangeMultiplier(2)
    ->Range(32, 2<<25)
    ->Complexity(benchmark::oN);
BENCHMARK(time_sparse_iteration_dense_vector)
    ->RangeMultiplier(2)
    ->Range(32, 2<<25)
    ->Complexity(benchmark::oN);