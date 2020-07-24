#include <UnitTest++/UnitTest++.h>
#include <libalgebra/vectors/vectors.h>

#include "time_and_details.h"

#include <string>
#include <map>
#include <sstream>
#include <iostream>


#define MAP_TYPE_PRINT() 

class Basis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef char KEY;


    typedef std::map<char, double> MAP;

    // Default constructor
    Basis() {}

    KEY begin() const 
    {
        return 'a';
    }

    KEY nextkey(const KEY& k) const
    {
        return k + 1;
    }

    KEY end() const
    {
        return 'z';
    }

    unsigned degree(const KEY &k) {
        return 1;
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<Basis*, KEY> &t
    ) {
        return os << t.second;
    }

};

typedef alg::vectors::dense_vector<Basis, 5> Vec;



SUITE(dense_vector_tests) {

    TEST(test_inserting_coordinate) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {};
        
        vec['a'] = 1.0;
        CHECK_EQUAL(1.0, vec['a']);   
    }

    TEST(text_formatting) {
        MAP_TYPE_PRINT();        
        TEST_DETAILS();
        std::basic_ostringstream<char> stream {};
        Vec vec {'a', 1.0};

        stream << vec;

        CHECK_EQUAL("{ 1(a) }", stream.str());
    }

    TEST(test_unary_minus) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        vec['b'] = 2.0;

        Vec nvec = -vec;
        CHECK_EQUAL(-1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_inplace_scalar_multiplication) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        double scalar = 2.0;

        vec *= scalar;
        CHECK_EQUAL(2.0, vec['a']);
    }    
    
    TEST(test_inplace_scalar_multiplication_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec neut {};
        double scalar = 2.0;

        neut *= scalar;

        CHECK_EQUAL(Vec {}, neut);
    }

    TEST(test_derived_scalar_multiplication) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        double scalar = 2.0;

        Vec nvec = vec*scalar;

        CHECK_EQUAL(2.0, nvec['a']);
    }

        
    TEST(test_derived_scalar_multiplication_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec neut {};
        double scalar = 2.0;

        Vec nvec = neut * scalar;

        CHECK_EQUAL(neut, nvec);
    }


    TEST(test_inplace_addition_same_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(3.0, vec1['a']);
    }

    TEST(test_inplace_addition_different_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(2.0, vec1['b']);
    }

    TEST(test_inplace_addition_interlaced_order) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'b', 1.0};
        vec1['d'] = 1.0;
        Vec vec2 {'a', 1.0};
        vec2['c'] = 1.0;

        Vec expected {'a', 1.0};
        expected['b'] = 1.0;
        expected['c'] = 1.0;
        expected['d'] = 1.0;

        vec1 += vec2;
        CHECK_EQUAL(expected, vec1);
    }

    TEST(test_inplace_addition_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        vec += neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_binary_addition_same_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        Vec nvec = vec1 + vec2;

        CHECK_EQUAL(3.0, nvec['a']);
    }

    TEST(test_derived_addition_different_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        Vec nvec = vec1 + vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(2.0, nvec['b']);
    }

    TEST(test_binary_addition_interlaced_order) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'b', 1.0};
        vec1['d'] = 1.0;
        Vec vec2 {'a', 1.0};
        vec2['c'] = 1.0;

        Vec expected {'a', 1.0};
        expected['b'] = 1.0;
        expected['c'] = 1.0;
        expected['d'] = 1.0;

        Vec nvec = vec1 + vec2;
        CHECK_EQUAL(expected, nvec);
    }
        
    TEST(test_derived_addition_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        Vec nvec = vec + neut;
        CHECK_EQUAL(1.0, nvec['a']);
    }

    TEST(test_inplace_subtraction_same_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(-1.0, vec1['a']);
    }

    TEST(test_inplace_subtract_same_keys_to_0) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        vec1 -= vec2;

        CHECK_EQUAL(0.0, vec1['a']);
    }

    TEST(test_inplace_subtraction_different_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(-2.0, vec1['b']);
    }

    TEST(test_inplace_subtraction_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        vec -= neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_derived_subtraction_same_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(-1.0, nvec['a']);
    }

    TEST(test_derived_subtraction_same_keys_to_0) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(0.0, nvec['a']);
    }

    TEST(test_derived_subtraction_different_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_derived_subtraction_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        Vec nvec = vec - neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_equality_equal_vectors) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        CHECK_EQUAL(vec1, vec2);
    }

    TEST(test_equality_not_equal_keys) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 1.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_not_equal_values) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_neutral_element) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();
        Vec neut1 {}, neut2 {};

        CHECK_EQUAL(neut1, neut2);
    }
    
    TEST(l1_norm_calculation) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();

        Vec vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;
        double expected = 6.0;

        CHECK_EQUAL(expected, vec.NormL1());
    }

    TEST(l1_norm_calculation_with_degree_1) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();

        Vec vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;

        // All elements in this basis have degree 1.
        double expected = 6.0;

        CHECK_EQUAL(expected, vec.NormL1(1));
    }

    TEST(l1_norm_calculation_with_degree_2) {
        MAP_TYPE_PRINT();
        TEST_DETAILS();

        Vec vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;

        // All elements in this basis have degree 1.
        double expected = 0.0;

        CHECK_EQUAL(expected, vec.NormL1(2));
    }




}


