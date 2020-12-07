//
// Created by sam on 01/12/2020.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include "time_and_details.h"
#include "alg_framework.h"
#include "helpers.h"

#include <random>
#include <unordered_set>

using alg::DEG;
using alg::DIMN;

template <DEG WIDTH, DEG DEPTH>
class InnovativePath
{
    typedef std::pair<int, int> coeff;


    std::vector<std::pair<size_t, coeff>> increments;

    template <typename S>
    std::vector<alg::lie<S, S, WIDTH, DEPTH>> make_increments() const
    {
        std::vector<alg::lie<S, S, WIDTH, DEPTH>> incr;
        incr.reserve(WIDTH);
        for (auto& i : increments) {
            S scalar {S{i.second.first} / S{i.second.second}};
            alg::LET letter {i.first + 1};
            incr.emplace_back(letter, scalar);
        }
        return incr;
    }

public:

    InnovativePath() : increments{}
    {
        std::mt19937 rng;
        rng.seed(12345);
        std::uniform_int_distribution<int> num_dist(-5, 5);
        std::uniform_int_distribution<int> den_dist(1, 8);

        increments.reserve(WIDTH);

        for (size_t i=0; i<WIDTH; ++i) {
            increments.emplace_back(i, coeff {num_dist(rng), den_dist(rng)});
        }
    }

    template <typename S>
    alg::free_tensor<S, S, WIDTH, DEPTH> signature() const
    {
        typedef alg::free_tensor<S, S, WIDTH, DEPTH> TENSOR;
        TENSOR sig(S(1));

        auto incr = make_increments<S>();
        alg::maps<S, S, WIDTH, DEPTH> maps;

        for (auto& inc : incr)
            sig *= exp(maps.l2t(inc));

        return sig;
    }

    template <typename S>
    alg::lie<S, S, WIDTH, DEPTH> log_signature() const
    {
        typedef alg::lie<S, S, WIDTH, DEPTH> LIE;
        alg::maps<S, S, WIDTH, DEPTH> maps;

        return maps.t2l(log(signature<S>()));
    }



};


template <typename S>
mpq_class scalar_error(mpq_class expected, S actual)
{
    return abs(mpq_class{actual} - expected);
}


SUITE(accuracy) {

    TEST(accuracy_signature_width_5_depth_5) {
        TEST_DETAILS();

        InnovativePath<5, 5> path;

        auto r_sig = path.signature<mpq_class>();
        auto f_sig = path.signature<float>();
        auto d_sig = path.signature<double>();

        auto r_it = r_sig.begin();
        auto f_it = f_sig.begin();
        auto d_it = d_sig.begin();

        auto test = [&] () { return r_it != r_sig.end() && f_it != f_sig.end()
                            && d_it != d_sig.end(); };

        for (; test(); ++r_it, ++f_it, ++d_it) {
            mpq_class f_err = scalar_error(r_it->second, f_it->second);
            mpq_class d_err = scalar_error(r_it->second, d_it->second);
            CHECK(f_err < 1e-5f);
            CHECK(d_err < 2e-10);
            CHECK_EQUAL(r_it->first, f_it->first);
            CHECK_EQUAL(r_it->first, d_it->first);
        }

    }



}