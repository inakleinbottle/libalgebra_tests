#include <benchmark/benchmark.h>

#include <libalgebra/vectors/vectors.h>
#include <random>
#include <iostream>

#ifndef USE_FLATMAP
#define MAP_T std::map<KEY, SCA>
#else
#include <boost/container/flat_map.hpp>
#define MAP_T boost::container::flat_map<KEY, SCA>
#endif


class Basis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef size_t KEY;
    typedef std::map<size_t, double> MAP;

    // Default constructor
    Basis() {}

    KEY begin() const 
    {
        return 0;
    }

    KEY nextkey(const KEY& k) const
    {
        return k + 1;
    }

    unsigned degree(const KEY &k) {
        return 1;
    }

    KEY end() const
    {
        return std::numeric_limits<KEY>::max();
    }

    inline static bool comp(const KEY& k1, const KEY& k2)
    {
        return k1 <= k2;
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<Basis*, KEY> &t
    ) {
        return os << t.second;
    }

};

template <size_t D>
using Vec = alg::vectors::dense_vector<Basis,D>;

template <size_t D>
Vec<D> generate_dense_vec(size_t num, size_t max_diff) {

    std::mt19937 key_rng;
    std::mt19937 sca_rng;

    std::uniform_real_distribution<double> reals (-2.0, 2.0);

    Vec<D> nvec;
    size_t crr = 0;
    size_t rn;
    for (size_t i=0; i<num; ++i) {
        //crr = crr + key_rng()%max_diff;
        //if (crr >= 25*num) throw "Bad index";
        try {
            nvec[i] = reals(sca_rng);
        } catch (alg::vectors::KeyNotFoundError) {
            std::cout << "missing key " << i << "\n";
            size_t k =0;
            for (auto u : nvec.get_keys()) {
                if (k == 5) break;
                std::cout << " " << u;
            }
            std::cout << "\n";
            throw;
        }
    }

    return nvec;
}


template<size_t D>
void time_addition_dense_vectors_t(benchmark::State &state)
{   
    Vec<D> vec1 = generate_dense_vec<D>(state.range(0), 25);
    Vec<D> vec2 = generate_dense_vec<D>(state.range(0), 25);

    for (auto _ : state) {
        benchmark::DoNotOptimize(vec1+vec2);
    }


}

#define CASE(i) case i : time_addition_dense_vectors_t<(i)>(state); break

static void time_addition_dense_vectors(benchmark::State &state)
{
    switch (state.range(0)) {
        CASE(2<<6);
        CASE(2<<7);
        CASE(2<<8);
        CASE(2<<9);
        CASE(2<<10);
        CASE(2<<11);
        CASE(2<<12);
        CASE(2<<13);
        CASE(2<<14);
        CASE(2<<15);
    }
}


BENCHMARK(time_addition_dense_vectors)->Range(2<<6, 2<<15);
