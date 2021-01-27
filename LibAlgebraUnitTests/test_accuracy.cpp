//
// Created by sam on 01/12/2020.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include "time_and_details.h"
#include "alg_framework.h"
#include "helpers.h"
#include "accuracy_test.h"

#include <random>
#include <unordered_set>





SUITE(accuracy) {

    using Framework = AccuracyTest<5, 5>;


    TEST_FIXTURE(Framework, accuracy_signature_width_5_depth_5) {
        TEST_DETAILS();

        std::mt19937 rng;
        uniform dist(-5, 5);
        std::vector<DIMN> seen { 1, 3, 2, 0, 4 };

        DIMN i = 0;

        auto filler = [&] (DIMN j) {
            if (j == seen[i])
                return Sca(1);
            return Sca(0);
        };

        std::vector<Lie> increments(5);
        for (; i<5; ++i)
            increments[i].fill_with(filler, 1);

        Maps fr_maps;

        TKey k;
        Tensor result (k);
        for (i=0; i<5; ++i)
            result *= exp(fr_maps.l2t(increments[i]));

        check(result);

    }



}