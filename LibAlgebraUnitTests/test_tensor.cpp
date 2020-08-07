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



SUITE(test_tensor_key_iterator) {

    typedef alg::free_tensor_basis<double, double, 2, 3> TBASIS;
    typedef typename TBASIS::KEY KEY;
    typedef alg::vectors::sparse_vector<TBASIS, std::map<KEY, double>> TENSOR;
    typedef typename TENSOR::iterator Iter;

    struct Fixture
    {
        TENSOR tensa, tensb, tensc;

        Fixture()
            : tensa {
                {KEY{}, 1.0}, {KEY{1}, 2.0}, {KEY{1, 1}, 3.0}
              },
              tensb {
                  {KEY{1}, 1.0}, {KEY{2}, 2.0}, {KEY{1, 1}, 3.0}, {KEY{1, 2}, 4.0}
              },
              tensc {
                  {KEY{}, 1.0}, {KEY{1,1}, 2.0}, {KEY{1,2}, 3.0}, {KEY{2,1}, 4.0}, {KEY{2, 2}, 5.0}
              }
        {}
    };

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_0) {
        TEST_DETAILS();
        Iter start = tensa.deg_begin(0);
        Iter end = tensa.deg_end(0);

        CHECK(start != end);
        CHECK(start == tensa.begin());
        CHECK(end != tensa.end());
        CHECK_EQUAL(KEY{}, start->first);
        CHECK_EQUAL(1.0, start->second);
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_0_not_present) {
        TEST_DETAILS();
        Iter start = tensb.deg_begin(0);
        Iter end = tensb.deg_end(0);

        Iter tend = tensb.end();
        CHECK(start == end);
        CHECK(start == tensb.begin());
        CHECK(end != tend);
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_1_not_start) {
        TEST_DETAILS();
        Iter start = tensa.deg_begin(1);
        Iter end = tensa.deg_end(1);

        CHECK(start != end);
        CHECK(start != tensa.begin());
        CHECK(start == ++tensa.begin());
        CHECK(end != tensa.end());
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_1_start) {
        TEST_DETAILS();
        Iter start = tensb.deg_begin(1);
        Iter end = tensb.deg_end(1);

        CHECK(start != end);
        CHECK(start == tensb.begin());
        CHECK(end != tensb.end());
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_1_missing) {
        TEST_DETAILS();
        Iter start = tensc.deg_begin(1);
        Iter end = tensc.deg_end(1);

        Iter tbegin = tensc.begin();
        CHECK(start == end);
        CHECK(start != tbegin);
        ++tbegin;
        CHECK(start == tbegin);
        CHECK(end != tensc.end());
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_2_a) {
        TEST_DETAILS();
        Iter start = tensa.deg_begin(2);
        Iter end = tensa.deg_end(2);

        CHECK(start != end);
        CHECK(start != tensa.begin());
        CHECK(start == ++(++tensa.begin()));
        CHECK(end == tensa.end());
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_2_b) {
        TEST_DETAILS();
        Iter start = tensb.deg_begin(2);
        Iter end = tensb.deg_end(2);

        CHECK(start != end);
        CHECK(start != tensb.begin());
        CHECK(start == ++(++tensb.begin()));
        CHECK(end == tensb.end());
    }

    TEST_FIXTURE(Fixture, test_degree_iterator_deg_2_c) {
        TEST_DETAILS();
        Iter start = tensc.deg_begin(2);
        Iter end = tensc.deg_end(2);

        CHECK(start != end);
        CHECK(start != tensc.begin());
        CHECK(start == ++tensc.begin());
        CHECK(end == tensc.end());
    }

    TEST_FIXTURE(Fixture, test_get_degree_0_then_degree_2) {
        TEST_DETAILS();

        Iter start = tensa.deg_begin(0);
        Iter end = tensa.deg_end(0);

        CHECK(start != end);
        CHECK(start == tensa.begin());
        CHECK(end != tensa.end());
        CHECK_EQUAL(KEY(), start->first);
        CHECK_EQUAL(1.0, start->second);

        start = tensa.deg_begin(2);
        

        CHECK(start != end);
        end = tensa.deg_end(2);

        CHECK(end == tensa.end());
        CHECK_EQUAL((KEY{1, 1}), start->first);
        CHECK_EQUAL(3.0, start->second);
    }

    TEST_FIXTURE(Fixture, test_iterator_keys_order_a) {
        TEST_DETAILS();

        Iter it = tensa.begin();
        Iter end = tensa.end();
        KEY k{};

        for (; it!=end; ++it) {
            CHECK(k == it->first || k < it->first);
            k = it->first;
        }
    }


    TEST_FIXTURE(Fixture, test_iterator_keys_order_b) {
        TEST_DETAILS();

        Iter it = tensb.begin();
        Iter end = tensb.end();
        KEY k{};

        for (; it!=end; ++it) {
            CHECK(k == it->first || k < it->first);
            k = it->first;
        }
    }


    TEST_FIXTURE(Fixture, test_iterator_keys_order_c) {
        TEST_DETAILS();

        Iter it = tensc.begin();
        Iter end = tensc.end();
        KEY k{};

        for (; it!=end; ++it) {
            CHECK(k == it->first || k < it->first);
            k = it->first;
        }
    }

}