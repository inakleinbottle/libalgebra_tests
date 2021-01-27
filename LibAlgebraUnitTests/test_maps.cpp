
#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>


#include <map>
#include <cstddef>
#include <string>
#include <cmath>

#include "time_and_details.h"
#include "alg_framework.h"
#include "helpers.h"





SUITE(maps_tests) {

    typedef alg_framework<4, 2, DPReal> framework;
    typedef typename framework::TENSOR TENSOR;
    typedef typename framework::LIE LIE;
    typedef typename TENSOR::KEY TKEY;

    TEST_FIXTURE(framework, test_expand_lie_key) {
        TEST_DETAILS();

        LET letter = 3; /// [1, 2]
        TKEY lhk, rhk;
        lhk.push_back(LET(1));
        lhk.push_back(LET(2));
        rhk.push_back(LET(2));
        rhk.push_back(LET(1));

        TENSOR expected {{lhk, 1.0}, {rhk, -1.0}};
        auto result = maps.expand(letter);
        CHECK_EQUAL(expected, result);

    }


    TEST_FIXTURE(framework, test_exp_log_roundtrip) {
        TEST_DETAILS();
        LET letter = 1;

        TENSOR exp_letter = maps.exp(letter);
        TENSOR rt_letter = log(exp_letter);

        TENSOR expected {letter, 1.0};

        CHECK_VEC_CLOSE(expected, rt_letter, 2.0e-15);
    }

    TEST_FIXTURE(framework, test_lie_to_tensor_single_letter) {
        TEST_DETAILS();
        LIE lie {1, 1.0};

        TENSOR expected {1, 1.0};

        CHECK_EQUAL(expected, maps.l2t(lie));
    }

    TEST_FIXTURE(framework, test_l2t_t2l_roundtrip) {
        TEST_DETAILS();

        LIE lie {{1, 1.0}, {2, 2.0}};

        CHECK_EQUAL(lie, maps.t2l(maps.l2t(lie)));

    }





}