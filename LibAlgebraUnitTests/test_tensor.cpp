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
#include "accuracy_test.h"

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


SUITE(tensor_basis) {

    typedef alg::free_tensor_basis<double, double, 2, 3> Basis;
    typedef typename Basis::KEY KEY;

    TEST(test_start_of_degree_0) {
        TEST_DETAILS();
        Basis b;
        CHECK_EQUAL(0, b.start_of_degree(0));
    }
        
    TEST(test_start_of_degree_1) {
        TEST_DETAILS();
        Basis b;
        CHECK_EQUAL(1, b.start_of_degree(1));
    }

    TEST(test_start_of_degree_2) {
        TEST_DETAILS();
        Basis b;
        CHECK_EQUAL(3, b.start_of_degree(2));
    }

    TEST(test_start_of_degree_3) {
        TEST_DETAILS();
        Basis b;
        CHECK_EQUAL(7, b.start_of_degree(3));
    }

    TEST(test_index_of_key_degree_0) {
        TEST_DETAILS();
        Basis b;
        KEY k {};

        CHECK_EQUAL(0, b.index_of_key(k));
    }

    TEST(test_index_of_key_degree_1) {
        TEST_DETAILS();
        Basis b;

        alg::DEG idx = 1;
        for (alg::LET i=1; i<=2; ++i) {
            KEY k {i};
            CHECK_EQUAL(idx, b.index_of_key(k));
            ++idx;
        }
    }

    TEST(test_index_of_key_degree_2) {
        TEST_DETAILS();
        Basis b;

        alg::DEG idx = 3;
        for (alg::LET i=1; i<=2; ++i) {
            for (alg::LET j=1; j<=2; ++j) {
                KEY k {i, j};
                CHECK_EQUAL(idx, b.index_of_key(k));
                ++idx;
            }

        }
    }

    TEST(test_index_of_key_degree_8) {
        TEST_DETAILS();
        Basis b;

        alg::DEG idx = 7;
        for (alg::LET i=1; i<=2; ++i) {
            for (alg::LET j=1; j<=2; ++j) {
                for (alg::LET k=1; k<=2; ++k) {
                    KEY key {i, j, k};
                    CHECK_EQUAL(idx, b.index_of_key(key));
                    ++idx;
                }
            }

        }
    }

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

SUITE(tensor_multiplication) {

    typedef alg_framework<5, 2, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef typename framework::SCA S;

    typedef alg::LET LET;

    TEST(test_multiplication_by_key_unit_unit_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_key(TKEY(), 1.0, TKEY(), 1.0, (typename TENSOR::VECT&)
        result);
        TENSOR expected {{TKEY(), 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_letter_unit_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_key(TKEY{1}, 1.0, TKEY(), 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_unit_letter_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_key(TKEY(), 1.0, TKEY{1}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_letter_letter_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_key(TKEY{1}, 1.0, TKEY{2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_unit_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_key(TKEY{1,1}, 1.0, TKEY{}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_unit_word_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_key(TKEY{}, 1.0, TKEY{1,1}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_letter_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(15)};
        op.apply_by_key(TKEY{1,1}, 1.0, TKEY{2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,1,2}, 1.0}};

                CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_letter_word_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(15)};
        op.apply_by_key(TKEY{1}, 1.0, TKEY{1,2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,1,2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_word_dense)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(15)};
        op.apply_by_key(TKEY{1,2}, 1.0, TKEY{1,2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,2,1,2}, 1.0}};

                CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_unit_unit_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY(), 1.0, TKEY(), 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY(), 1.0}};

        CHECK_EQUAL(expected, result);
    }


    TEST(test_multiplication_by_key_letter_unit_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1}, 1.0, TKEY(), 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1}, 1.0}};

                CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_unit_letter_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY(), 1.0, TKEY{1}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_letter_letter_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1}, 1.0, TKEY{2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_unit_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1,1}, 1.0, TKEY{}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_unit_word_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{}, 1.0, TKEY{1,1}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_letter_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1,1}, 1.0, TKEY{2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,1,2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_letter_word_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1}, 1.0, TKEY{1,2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,1,2}, 1.0}};

                CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_key_word_word_sparse)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{};
        op.apply_by_key(TKEY{1,2}, 1.0, TKEY{1,2}, 1.0, (typename TENSOR::VECT&)
                result);
        TENSOR expected {{TKEY{1,2,1,2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

 /*   TEST(test_multiplication_by_index_unit_unit)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_index(0, 1.0, 0, 0, 1.0, (typename TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY(), 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_unit_letter)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_index(0, 1.0, 1, 1, 1.0, (typename TENSOR::VECT&) result);
        TENSOR expected {{alg::LET{1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_letter_unit)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(3)};
        op.apply_by_index(1, 1.0, 0, 0, 1.0, (typename TENSOR::VECT&) result);
        TENSOR expected {{alg::LET{1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }


    TEST(test_multiplication_by_index_letter_letter)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_index(1, 1.0, 1, 1, 1.0, (typename TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_word_unit)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_index(3 *//*(1,1)*//*, 1.0, 0, 0, 1.0, (typename
        TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_unit_word)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(7)};
        op.apply_by_index(0, 1.0, 2, 3 *//*(1,1)*//*, 1.0, (typename
        TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1}, 1.0}};

                CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_word_letter)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(15)};
        op.apply_by_index(3 *//*(1,1)*//*, 1.0, 1, 1, 1.0, (typename
        TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_letter_word)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(15)};
        op.apply_by_index(1, 1.0, 2, 3*//*(1,1)*//*, 1.0, (typename
        TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1, 1}, 1.0}};

        CHECK_EQUAL(expected, result);
    }

    TEST(test_multiplication_by_index_word_word)
    {
        TEST_DETAILS();

        typedef typename TENSOR::ALG::scalar_passthrough SPT;
        typename TENSOR::BASIS::MultiplyBasisOperator<SPT> op{SPT()};

        TENSOR result{TENSOR::VECT::create_with_dimension(31)};
        op.apply_by_index(3 *//*(1,1)*//*, 1.0, 2 , 4*//*(1,2)*//*, 1.0, (typename
        TENSOR::VECT&) result);
        TENSOR expected {{typename TENSOR::KEY {1, 1, 1, 2}, 1.0}};

        CHECK_EQUAL(expected, result);
    }*/
}

SUITE(tensor_multiplication_2_5_double) {
    typedef alg_framework<5, 2, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef alg::LET LET;
    typedef typename framework::SCA S;


    TEST(test_multiplication_dense_low_dimension)
    {
        TEST_DETAILS();

        TENSOR t1 {
                {TKEY {}, 1.0},
                {TKEY {1}, 2.0},
                {TKEY {2}, 3.0}
        };

        TENSOR t2 {
                {TKEY {}, -1.0},
                {TKEY {1}, -2.0},
                {TKEY {2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        TENSOR expected {
                {TKEY {}, -1.0},
                {TKEY {1}, -4.0},
                {TKEY {2}, -6.0},
                {TKEY {1, 1}, -4.0},
                {TKEY {1, 2}, -6.0},
                {TKEY {2, 1}, -6.0},
                {TKEY {2, 2}, -9.0}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

    }

    TEST(test_multiplication_dense_and_sparse)
    {
        TEST_DETAILS();

        TENSOR t1 {
                {TKEY {}, 1.0},
                {TKEY {1}, 2.0},
                {TKEY {2}, 3.0}
        };

        TENSOR t2 {
                {TKEY {}, -1.0},
                {TKEY {1}, -2.0},
                {TKEY {2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        t1[TKEY {1,2}] = 0.5;

        TENSOR expected {
                {TKEY {}, -1.0},
                {TKEY {1}, -4.0},
                {TKEY {2}, -6.0},
                {TKEY {1, 1}, -4.0},
                {TKEY {1, 2}, -6.5},
                {TKEY {2, 1}, -6.0},
                {TKEY {2, 2}, -9.0},
                {TKEY {1, 2, 1}, -1.0},
                {TKEY {1, 2, 2}, -1.5}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

    }

}


SUITE(tensor_multiplication_5_5_double) {
    typedef alg_framework<5, 5, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef alg::LET LET;
    typedef typename framework::SCA S;


    TEST (test_multiplication_dense_low_dimension) {
        TEST_DETAILS();

        TENSOR t1{
                {TKEY{},  1.0},
                {TKEY{1}, 2.0},
                {TKEY{2}, 3.0}
        };

        TENSOR t2{
                {TKEY{},  -1.0},
                {TKEY{1}, -2.0},
                {TKEY{2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        TENSOR expected{
                {TKEY{},     -1.0},
                {TKEY{1},    -4.0},
                {TKEY{2},    -6.0},
                {TKEY{1, 1}, -4.0},
                {TKEY{1, 2}, -6.0},
                {TKEY{2, 1}, -6.0},
                {TKEY{2, 2}, -9.0}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

    }

    TEST (test_multiplication_dense_and_sparse) {
        TEST_DETAILS();

        TENSOR t1{
                {TKEY{},  1.0},
                {TKEY{1}, 2.0},
                {TKEY{2}, 3.0}
        };

        TENSOR t2{
                {TKEY{},  -1.0},
                {TKEY{1}, -2.0},
                {TKEY{2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        t1[TKEY{1, 2}] = 0.5;

        TENSOR expected{
                {TKEY{},        -1.0},
                {TKEY{1},       -4.0},
                {TKEY{2},       -6.0},
                {TKEY{1, 1},    -4.0},
                {TKEY{1, 2},    -6.5},
                {TKEY{2, 1},    -6.0},
                {TKEY{2, 2},    -9.0},
                {TKEY{1, 2, 1}, -1.0},
                {TKEY{1, 2, 2}, -1.5}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

    }

}


SUITE(tensor_multiplication_8_5_double) {
    typedef alg_framework<5, 8, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::TENSOR::KEY TKEY;
    typedef alg::LET LET;
    typedef typename framework::SCA S;


    TEST(test_multiplication_dense_low_dimension)
    {
        TEST_DETAILS();

        TENSOR t1 {
                {TKEY {}, 1.0},
                {TKEY {1}, 2.0},
                {TKEY {2}, 3.0}
        };

        TENSOR t2 {
                {TKEY {}, -1.0},
                {TKEY {1}, -2.0},
                {TKEY {2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        TENSOR expected {
                {TKEY {}, -1.0},
                {TKEY {1}, -4.0},
                {TKEY {2}, -6.0},
                {TKEY {1, 1}, -4.0},
                {TKEY {1, 2}, -6.0},
                {TKEY {2, 1}, -6.0},
                {TKEY {2, 2}, -9.0}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

    }

    TEST(test_multiplication_dense_and_sparse)
    {
        TEST_DETAILS();

        TENSOR t1 {
                {TKEY {}, 1.0},
                {TKEY {1}, 2.0},
                {TKEY {2}, 3.0}
        };

        TENSOR t2 {
                {TKEY {}, -1.0},
                {TKEY {1}, -2.0},
                {TKEY {2}, -3.0}
        };
        t1.maybe_resize(1);
        t2.maybe_resize(1);

        t1[TKEY {1,2}] = 0.5;

        TENSOR expected {
                {TKEY {}, -1.0},
                {TKEY {1}, -4.0},
                {TKEY {2}, -6.0},
                {TKEY {1, 1}, -4.0},
                {TKEY {1, 2}, -6.5},
                {TKEY {2, 1}, -6.0},
                {TKEY {2, 2}, -9.0},
                {TKEY {1, 2, 1}, -1.0},
                {TKEY {1, 2, 2}, -1.5}
        };

        TENSOR t3 = t1 * t2;

                CHECK_EQUAL(expected, t3);

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


SUITE(test_dense_tensor) {

    typedef alg::free_tensor_basis<double, double, 2, 3> TBASIS;
    typedef typename TBASIS::KEY KEY;
    typedef alg::vectors::dense_vector<TBASIS> TENSOR;


    TEST(test_basis_key_indices) {
        TEST_DETAILS();
        TENSOR tens (15);
        TBASIS b;
        constexpr size_t dimension = 1 + 2 + 2*2 + 2*2*2;

        double c = 0.0;
        for (KEY k = b.begin(); !(k == b.end()); k=b.nextkey(k)) {
            tens[k] = c;
            c += 1.0;
        }

        c = 0.0;
        bool skip = true;
        auto coeffs = tens.get_coeffs();
        for (size_t i=0; i < dimension; ++i) {
            CHECK_EQUAL(c, coeffs[i]);
//            if (skip) {
//                skip = false;
//                continue;
//            }
            c += 1.0;
        }
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

#define _TEST(name) TEST_FIXTURE(Fixture, name)

SUITE (tensor_ops_accuracy_tensor22) {

    using Fixture = AccuracyTest<2, 8>;

    _TEST(random_tensor_mul_max_deg) {
        TEST_DETAILS();

        Tensor a, b;
        std::mt19937 rng;
        uniform dist(-5, 5);

        auto fill = [&](DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        auto result = a * b;

        check(result);
    }

    _TEST(random_tensor_mul_max_deg_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        std::mt19937 rng;
        uniform dist(-5, 5);

        auto fill = [&](DIMN i) {
            if (i == 0) return Sca(1);
            return dist(rng);
        };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        auto result = a * b;

        check(result);
    }

    _TEST(random_tensor_imul_max_deg) {
        TEST_DETAILS();

        Tensor a, b;
        std::mt19937 rng;
        uniform dist(-5, 5);

        auto fill = [&](DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;

        check(a);
    }

    _TEST(random_tensor_imul_max_deg_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        std::mt19937 rng;
        uniform dist(-5, 5);

        auto fill = [&](DIMN i) {
            if (i == 0) return Sca(1);
            return dist(rng);
        };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;

        check(a);
    }
}


SUITE(tensor_ops_accuracy_tensor28) {

    using Fixture = AccuracyTest<2, 8>;

    _TEST(random_tensor_mul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_mul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_imul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 3);
        b.fill_with(fill, 3);

        a *= b;
        check(a);
    }

    _TEST(random_tensor_imul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) {
            if (i == 0) return Sca(1);
            return dist(rng);
        };
        a.fill_with(fill, 3);
        b.fill_with(fill, 3);

        a *= b;
        check(a);
    }

    _TEST(complete_random_mul_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        a.fill_with([&] (DIMN i) { return dist(rng); }, 2);
        b.fill_with([&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); }, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_exp_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_zero_elt) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_non_zero) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_log_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

    _TEST(random_tensor_inv_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

    _TEST(random_tensor_mul_with_sparse) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };

        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a.insert_sparse(fill, 3, 5, rng, dist);
        b.insert_sparse(fill, 3, 5, rng, dist);

        auto result = a * b;

        check(result);
    }

    _TEST(random_tensor_imul_with_sparse) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };

        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a.insert_sparse(fill, 3, 5, rng, dist);
        b.insert_sparse(fill, 3, 5, rng, dist);

        a *= b;

        check(a);
    }

    _TEST(random_tensor_mul_unit_1_with_sparse) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };

        a.fill_with(fill, 2);
        b.fill_with([&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); }, 2);

        a.insert_sparse(fill, 3, 5, rng, dist);
        b.insert_sparse(fill, 3, 5, rng, dist);

        auto result = a * b;

        check(a);
    }


    _TEST(random_tensor_exp_with_sparse) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };

        a.fill_with(fill, 2);

        a.insert_sparse(fill, 3, 5, rng, dist);

        auto result = exp(a);

        check(result, 2.0e-3, 2.0e-11);
    }

    _TEST(random_tensor_exp_deg_1_sparse) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };

        a.insert_sparse(fill, 1, 1, rng, dist, 2);

        auto result = exp(a);

        check(result, 2.0e-3, 2.0e-11);
    }



}


SUITE (tensor_ops_accuracy_tensor32) {

    using Fixture = AccuracyTest<3, 2>;

    _TEST(random_tensor_mul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_mul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_imul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;
        check(a);
    }

    _TEST(random_tensor_imul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) {
            if (i == 0) return Sca(1);
            return dist(rng);
        };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;
        check(a);
    }

    _TEST(complete_random_mul_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        a.fill_with([&] (DIMN i) { return dist(rng); }, 2);
        b.fill_with([&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); }, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_exp_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_zero_elt) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_non_zero) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_log_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

    _TEST(random_tensor_inv_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

}


SUITE (tensor_ops_accuracy_tensor82) {

    using Fixture = AccuracyTest<8, 2>;

    _TEST(random_tensor_mul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_mul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_imul_accuracy) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;
        check(a);
    }

    _TEST(random_tensor_imul_accuracy_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) {
            if (i == 0) return Sca(1);
            return dist(rng);
        };
        a.fill_with(fill, 2);
        b.fill_with(fill, 2);

        a *= b;
        check(a);
    }

    _TEST(complete_random_mul_unit_1) {
        TEST_DETAILS();

        Tensor a, b;
        uniform dist(-5, 5);
        std::mt19937 rng;

        a.fill_with([&] (DIMN i) { return dist(rng); }, 2);
        b.fill_with([&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); }, 2);

        Tensor result = a * b;

        check(result);
    }

    _TEST(random_tensor_exp_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_zero_elt) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_exp_accuracy_deg1_non_zero) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { return dist(rng); };
        a.fill_with(fill, 1);

        Tensor b = exp(a);

        check(b);
    }

    _TEST(random_tensor_log_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(0); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

    _TEST(random_tensor_inv_accuracy) {
        TEST_DETAILS();

        Tensor a;
        uniform dist(-5, 5);
        std::mt19937 rng;

        auto fill = [&] (DIMN i) { if (i==0) return Sca(1); return dist(rng); };
        a.fill_with(fill, 2);

        Tensor b = log(a);

        check(b);
    }

}


#undef _TEST