#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include <addons/gmpwrapper.h>

#include <map>
#include <cstddef>
#include <string>
#include <cmath>

#include "time_and_details.h"
#include "alg_framework.h"
#include "helpers.h"


template<class T, class S>
T exp_to_depth(T x, size_t depth, S one) {
    T result (one);
    T xn(one);
    size_t fact = 1;
    for (size_t i=1; i<=depth; ++i) {
        xn *= x;
        fact *= i;
        result += (xn / fact);
    }
    return result;
}


SUITE(free_tensor_tests_double) {

    typedef alg_framework<5, 2, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef typename framework::SCA S;
    typedef alg::LET LET;

    TEST(tensor_exponential_zero) {
        TEST_DETAILS();
        TENSOR ten(S(0));
        TENSOR expected {TKEY()};

        CHECK_EQUAL(expected, exp(ten));
    }

    TEST(test_exponential_tensor_unit) {
        TEST_DETAILS();

        TENSOR ten {TKEY()};
        TENSOR expected(exp_to_depth(S(1), 5, S(1)));
        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_single_letter) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten {{letter, S(1)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_multiple_letter) {
        TEST_DETAILS();

        TENSOR ten {{LET(1), S(1)}, {LET(2), S(1)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_single_letter_with_coeff) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten {{letter, S(2)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);        
    }

    TEST(test_exponential_multiple_letter_with_coeffs) {
        TEST_DETAILS();

        TENSOR ten {{LET(1), S(1)}, {LET(2), S(2)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_log_tensor_unit) {
        TEST_DETAILS();

        TENSOR tunit {TKEY()};
        TENSOR zero {S(0)};

        CHECK_EQUAL(zero, log(tunit));
    }

    TEST(test_log_tensor_with_no_explicit_unit) {
        TEST_DETAILS();

        TENSOR tunit {TKEY()};
        TENSOR no_unit {LET(1), S(1)};
        TENSOR ten = tunit + no_unit;

        CHECK_EQUAL(log(ten), log(no_unit));
    }

    TEST(test_log_exp_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten {LET(1), S(1)};
        
        CHECK_VEC_CLOSE(ten, log(exp(ten)), 2.0e-15);
    }

    TEST(test_exp_log_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten {LET(1), S(1)};
        TENSOR expected = TENSOR{TKEY()} + ten;

        CHECK_VEC_CLOSE(expected, exp(log(ten)), 2.0e-15);
    }

}



SUITE(test_free_tensor_rational) {

    typedef alg_framework<5, 2, Rational> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef alg::LET LET;
    typedef typename framework::SCA S;

    TEST(tensor_exponential_zero) {
        TEST_DETAILS();
        TENSOR ten(S(0));
        TENSOR expected {TKEY()};

        CHECK_EQUAL(expected, exp(ten));
    }

    TEST(test_exponential_tensor_unit) {
        TEST_DETAILS();

        TENSOR ten {TKEY()};
        TENSOR expected(exp_to_depth(S(1), 5, S(1)));
        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_single_letter) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten {{letter, S(1)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_multiple_letter) {
        TEST_DETAILS();

        TENSOR ten {{LET(1), S(1)}, {LET(2), S(1)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_exponential_single_letter_with_coeff) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR ten {{letter, S(2)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);        
    }

    TEST(test_exponential_multiple_letter_with_coeffs) {
        TEST_DETAILS();

        TENSOR ten {{LET(1), S(1)}, {LET(2), S(2)}};
        TENSOR expected = exp_to_depth(ten, 5, S(1));

        CHECK_VEC_CLOSE(expected, exp(ten), 2.0e-15);
    }

    TEST(test_log_tensor_unit) {
        TEST_DETAILS();

        TENSOR tunit {TKEY()};
        TENSOR zero {S(0)};

        CHECK_EQUAL(zero, log(tunit));
    }

    TEST(test_log_tensor_with_no_explicit_unit) {
        TEST_DETAILS();

        TENSOR tunit {TKEY()};
        TENSOR no_unit {LET(1), S(1)};
        TENSOR ten = tunit + no_unit;

        CHECK_EQUAL(log(ten), log(no_unit));
    }

    TEST(test_log_exp_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten {LET(1), S(1)};
        
        CHECK_VEC_CLOSE(ten, log(exp(ten)), 2.0e-15);
    }

    TEST(test_exp_log_round_trip_single_letter) {
        TEST_DETAILS();

        TENSOR ten {LET(1), S(1)};
        TENSOR expected = TENSOR{TKEY()} + ten;

        CHECK_VEC_CLOSE(expected, exp(log(ten)), 2.0e-15);
    }


}