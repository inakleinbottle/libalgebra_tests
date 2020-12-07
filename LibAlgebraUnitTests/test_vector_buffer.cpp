//
// Created by sam on 23/11/2020.
//

#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>


#include "time_and_details.h"

#include <map>
#include <utility>
#include <iostream>


#define WIDTH 5
#define DEPTH 8


typedef alg::DEG DEG;
typedef alg::DIMN DIMN;


SUITE(vector_buffer) {

    using Basis = alg::free_tensor_basis<float, float, WIDTH, DEPTH>;
    using VectorBuffer = alg::vectors::vector_buffer<
            Basis,
            typename Basis::SCALAR,
            std::allocator<typename Basis::SCALAR>,
            32 / sizeof(float)
            >;


    TEST(test_resize_deg_boundary_0)
    {
        TEST_DETAILS();
        VectorBuffer b;
        b.resize(1);

        CHECK_EQUAL(8, b.raw_size());
    }

    TEST(test_resize_deg_boundary_1)
    {
        TEST_DETAILS();
        VectorBuffer b;
        b.resize(1 + WIDTH);

        CHECK_EQUAL(16, b.raw_size());
    }

    TEST(test_resize_deg_boundary_2)
    {
        TEST_DETAILS();
        VectorBuffer b;
        b.resize(1 + WIDTH + WIDTH * WIDTH);

        CHECK_EQUAL(48, b.raw_size());
    }

    TEST(test_get_by_index_0_no_hint)
    {
        TEST_DETAILS();
        VectorBuffer b {
            {1.0f, 0.0f, 0.0f, 0.0f,
             0.0f, 0.0f, 0.0f, 0.0f},
        };

        CHECK_EQUAL(1.0f, b.get(0));
    }

    TEST(test_get_by_index_1_no_hint) {
        TEST_DETAILS();
        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                 6.0f, 0.0f, 0.0f, 0.0f},
        };

        CHECK_EQUAL(2.0f, b.get(1));
        CHECK_EQUAL(3.0f, b.get(2));
        CHECK_EQUAL(4.0f, b.get(3));
        CHECK_EQUAL(5.0f, b.get(4));
        CHECK_EQUAL(6.0f, b.get(5));

    }

    TEST(test_get_by_index_2_no_hint) {
        TEST_DETAILS();
        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                 6.0f, 0.0f, 0.0f, 0.0f},
                {7.0f, 8.0f, 9.0f, 10.0f,
                 11.0f, 12.0f, 13.0f, 14.0f},
                {15.0f, 16.0f, 17.0f, 18.0f,
                 19.0f, 20.0f, 21.0f, 22.0f},
                {23.0f, 24.0f, 25.0f, 26.0f,
                 27.0f, 28.0f, 29.0f, 30.0f},
                {31.0f, 32.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
        };

        for (DIMN i=0; i<25; ++i)
            CHECK_EQUAL(7.0f + i, b.get(6+i));
    }


    TEST(test_get_by_index_0_hint)
    {
        TEST_DETAILS();
        VectorBuffer b {
                {1.0f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
        };

                CHECK_EQUAL(1.0f, b.get(0, 0));
    }

    TEST(test_get_by_index_1_hint) {
        TEST_DETAILS();
        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                        6.0f, 0.0f, 0.0f, 0.0f},
        };

                CHECK_EQUAL(2.0f, b.get(1, 1));
                CHECK_EQUAL(3.0f, b.get(2, 1));
                CHECK_EQUAL(4.0f, b.get(3, 1));
                CHECK_EQUAL(5.0f, b.get(4, 1));
                CHECK_EQUAL(6.0f, b.get(5, 1));

    }

    TEST(test_get_by_index_2_hint) {
        TEST_DETAILS();
        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                        6.0f, 0.0f, 0.0f, 0.0f},
                {7.0f, 8.0f, 9.0f, 10.0f,
                        11.0f, 12.0f, 13.0f, 14.0f},
                {15.0f, 16.0f, 17.0f, 18.0f,
                        19.0f, 20.0f, 21.0f, 22.0f},
                {23.0f, 24.0f, 25.0f, 26.0f,
                        27.0f, 28.0f, 29.0f, 30.0f},
                {31.0f, 32.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
        };

        for (DIMN i=0; i<25; ++i)
            CHECK_EQUAL(7.0f + i, b.get(6+i, 2));
    }

    TEST(test_iterator_deg_0) {
        TEST_DETAILS();

        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
        };

        std::vector<float> expected {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

        auto cit = expected.cbegin();

        for (auto item : b)
            CHECK_EQUAL(*(cit++), item);

    }

    TEST(test_iterator_deg_1) {
        TEST_DETAILS();

        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
        };

        std::vector<float> expected {
            1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, // deg 0
            2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 0.0f, 0.0f, 0.0f, // deg 1
        };

        auto cit = expected.cbegin();

        for (auto item : b)
            CHECK_EQUAL(*(cit++), item);
    }


    TEST(test_deg_start_0) {
        TEST_DETAILS();

        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                 6.0f, 0.0f, 0.0f, 0.0f},
                {7.0f, 8.0f, 9.0f, 10.0f,
                 11.0f, 12.0f, 13.0f, 14.0f},
                {15.0f, 16.0f, 17.0f, 18.0f,
                 19.0f, 20.0f, 21.0f, 22.0f},
                {23.0f, 24.0f, 25.0f, 26.0f,
                 27.0f, 28.0f, 29.0f, 30.0f},
                {31.0f, 32.0f, 0.0f, 0.0f,
                 0.0f, 0.0f, 0.0f, 0.0f},
        };

        CHECK_EQUAL(1.0f, *b.deg_start(0));
    }

    TEST(test_deg_start_1) {
        TEST_DETAILS();

        VectorBuffer b{
                {1.0f, 0.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
                {2.0f, 3.0f, 4.0f, 5.0f,
                        6.0f, 0.0f, 0.0f, 0.0f},
                {7.0f, 8.0f, 9.0f, 10.0f,
                        11.0f, 12.0f, 13.0f, 14.0f},
                {15.0f, 16.0f, 17.0f, 18.0f,
                        19.0f, 20.0f, 21.0f, 22.0f},
                {23.0f, 24.0f, 25.0f, 26.0f,
                        27.0f, 28.0f, 29.0f, 30.0f},
                {31.0f, 32.0f, 0.0f, 0.0f,
                        0.0f, 0.0f, 0.0f, 0.0f},
        };

                CHECK_EQUAL(2.0f, *b.deg_start(1));
    }

    TEST(test_deg_start_2) {
        TEST_DETAILS();

        VectorBuffer b{
                {1.0f,  0.0f,  0.0f,  0.0f,
                        0.0f,  0.0f,  0.0f,  0.0f},
                {2.0f,  3.0f,  4.0f,  5.0f,
                        6.0f,  0.0f,  0.0f,  0.0f},
                {7.0f,  8.0f,  9.0f,  10.0f,
                        11.0f, 12.0f, 13.0f, 14.0f},
                {15.0f, 16.0f, 17.0f, 18.0f,
                        19.0f, 20.0f, 21.0f, 22.0f},
                {23.0f, 24.0f, 25.0f, 26.0f,
                        27.0f, 28.0f, 29.0f, 30.0f},
                {31.0f, 32.0f, 0.0f,  0.0f,
                        0.0f,  0.0f,  0.0f,  0.0f},
        };

        CHECK_EQUAL(7.0f, *b.deg_start(2));
    }

    TEST(test_degree_target_0)
    {
        VectorBuffer b;
        CHECK_EQUAL(1, b.deg_target(0));
    }

    TEST(test_degree_target_1)
    {
        VectorBuffer b;
        CHECK_EQUAL(5, b.deg_target(1));
    }

    TEST(test_degree_target_2)
    {
        VectorBuffer b;
        CHECK_EQUAL(25, b.deg_target(2));
    }
}