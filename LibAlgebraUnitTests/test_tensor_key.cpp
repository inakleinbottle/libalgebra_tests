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
#include <libalgebra/constpower.h>

SUITE(tensor_basis_key) {

    typedef alg::_tensor_basis<2, 5> KEY;

    TEST(test_initializer_list_single_key) {
        TEST_DETAILS();

        KEY key {{1}};
        KEY expected{1};

        CHECK_EQUAL(expected, key);
    }

    TEST(test_key_less_empty_key_degree_1) {
        TEST_DETAILS();

        KEY key1 {};
        KEY key2 {2};

        CHECK(key1 < key2);
    }

    TEST(test_key_less_same_degree_1) {
        TEST_DETAILS();

        KEY key1 {1};
        KEY key2 {2};

        CHECK(key1 < key2);
    }

    TEST(test_key_less_different_degree) {
        TEST_DETAILS();

        KEY key1 {1};
        KEY key2 {2, 2};

        CHECK(key1 < key2);
    }

    TEST(test_degree_2_less_ordering) {
        TEST_DETAILS();

        KEY k1 {1,1}, k2{1,2}, k3{2, 1}, k4{2, 2};

        CHECK(k1 < k2);
        CHECK(k2 < k3);
        CHECK(k3 < k4);
    }

    TEST(test_degree_3_degree_2_less) {
        TEST_DETAILS();

        KEY k1 {2, 2}, k2 {1, 1, 1};

        CHECK(k1 < k2);
    }

    TEST(test_firstletter_and_rparent) {
        TEST_DETAILS();

        KEY k1 {1, 2, 1};
        KEY k = k1.FirstLetter();
        KEY p = k1.rparent();

        CHECK_EQUAL(k1, k * p);
    }

    TEST(test_index_from_key_roundtrip) {
        TEST_DETAILS();

        alg::free_tensor_basis<double, double, 2, 5> basis;

        size_t basis_size = (alg::ConstPower<2, 6>::ans - 1);

        for (size_t i=0; i<basis_size; ++i) {
            KEY k = basis.key_of_index(i);
            CHECK_EQUAL(i, basis.index_of_key(k));
        }
    }

}