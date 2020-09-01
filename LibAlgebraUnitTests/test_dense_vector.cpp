#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include "time_and_details.h"

#include <string>
#include <map>
#include <sstream>
#include <iostream>

typedef alg::DEG DEG;
typedef alg::DIMN DIMN;

class DVTBasis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef char KEY;
    static const DEG MAX_DEGREE = 1;

    typedef std::map<char, double> MAP;

    // Default constructor
    DVTBasis() {}

    KEY begin() const 
    {
        return 'a';
    }

    KEY nextkey(const KEY& k) const
    {
        return KEY(k + 1);
    }

    KEY end() const
    {
        return 'z';
    }

    unsigned degree(const KEY &k) {
        return 1;
    }

    static constexpr DIMN max_dimension()
    {
        return 26;
    }

    inline static DIMN index_of_key(const KEY& k)
    {
        return static_cast<DIMN>(k - 'a');
    }

    inline static KEY key_of_index(const DIMN& idx)
    {
        return static_cast<KEY>(idx + 'a');
    }

    static DIMN start_of_degree(const DEG& d)
    {
        if (d <= 1)
            return 0;
        return 26;
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<DVTBasis*, KEY> &t
    ) {
        return os << t.second;
    }

};

typedef alg::vectors::dense_vector<DVTBasis> Vec;



SUITE(dense_vector_tests) {

    TEST(test_vector_indexing) {
        TEST_DETAILS();
        Vec vec {{'a', 1.0}, {'b', 2.0}, {'e', 5.0}};

        CHECK_EQUAL(1.0, vec['a']);
        CHECK_EQUAL(2.0, vec['b']);
        CHECK_EQUAL(0.0, vec['c']);
        CHECK_EQUAL(0.0, vec['d']);
        CHECK_EQUAL(5.0, vec['e']);
    }
    
    
    TEST(test_inserting_coordinate) {
        TEST_DETAILS();
        Vec vec(5);
        
        vec['a'] = 1.0;
        CHECK_EQUAL(1.0, vec['a']);   
    }

    TEST(text_formatting) {        
        TEST_DETAILS();
        std::basic_ostringstream<char> stream {};
        Vec vec {'a', 1.0};

        stream << vec;

        CHECK_EQUAL("{ 1(a) }", stream.str());
    }

    TEST(test_unary_minus) {
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        vec['b'] = 2.0;

        Vec nvec = -vec;
        CHECK_EQUAL(-1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_inplace_scalar_multiplication) {
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        double scalar = 2.0;

        vec *= scalar;
        CHECK_EQUAL(2.0, vec['a']);
    }    
    
    TEST(test_inplace_scalar_multiplication_neutral_element) {
        TEST_DETAILS();
        Vec neut {};
        double scalar = 2.0;

        neut *= scalar;

        CHECK_EQUAL(Vec {}, neut);
    }

    TEST(test_derived_scalar_multiplication) {
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        double scalar = 2.0;

        Vec nvec = vec*scalar;

        CHECK_EQUAL(2.0, nvec['a']);
    }

        
    TEST(test_derived_scalar_multiplication_neutral_element) {
        TEST_DETAILS();
        Vec neut {};
        double scalar = 2.0;

        Vec nvec = neut * scalar;

        CHECK_EQUAL(neut, nvec);
    }

    TEST(test_inplace_rational_division) {
        TEST_DETAILS();
        Vec vec = {{'a', 1.0}, {'b', 1.0}, {'c', 1.0}};
        double rational = 2.0;

        vec /= rational;
        Vec expected = {{'a', 0.5}, {'b', 0.5}, {'c', 0.5}};
        CHECK_EQUAL(expected, vec);
    }

    TEST(test_binary_rational_division) {
        TEST_DETAILS();
        Vec vec {{'a', 1.0}, {'b', 1.0}, {'c', 1.0}};
        double rational = 2.0;

        Vec expected = {{'a', 0.5}, {'b', 0.5}, {'c', 0.5}};
        CHECK_EQUAL(expected, vec / rational);
    }


    TEST(test_inplace_addition_same_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(3.0, vec1['a']);
    }

    TEST(test_inplace_addition_different_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(2.0, vec1['b']);
    }

    TEST(test_inplace_addition_interlaced_order) {
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
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        vec += neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_binary_addition_same_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        Vec nvec = vec1 + vec2;

        CHECK_EQUAL(3.0, nvec['a']);
    }

    TEST(test_derived_addition_different_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        Vec nvec = vec1 + vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(2.0, nvec['b']);
    }

    TEST(test_binary_addition_interlaced_order) {
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
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        Vec nvec = vec + neut;
        CHECK_EQUAL(1.0, nvec['a']);
    }

    TEST(test_inplace_subtraction_same_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(-1.0, vec1['a']);
    }

    TEST(test_inplace_subtract_same_keys_to_0) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        vec1 -= vec2;

        CHECK_EQUAL(0.0, vec1['a']);
    }

    TEST(test_inplace_subtraction_different_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(-2.0, vec1['b']);
    }

    TEST(test_inplace_subtraction_neutral_element) {
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        vec -= neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_derived_subtraction_same_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(-1.0, nvec['a']);
    }

    TEST(test_derived_subtraction_same_keys_to_0) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(0.0, nvec['a']);
    }

    TEST(test_derived_subtraction_different_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 2.0};

        Vec nvec = vec1 - vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_derived_subtraction_neutral_element) {
        TEST_DETAILS();
        Vec vec {'a', 1.0};
        Vec neut {};

        Vec nvec = vec - neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_equality_equal_vectors) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        CHECK_EQUAL(vec1, vec2);
    }

    TEST(test_equality_not_equal_keys) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 1.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_not_equal_values) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_neutral_element) {
        TEST_DETAILS();
        Vec neut1 {}, neut2 {};

        CHECK_EQUAL(neut1, neut2);
    }
    
    TEST(l1_norm_calculation) {
        TEST_DETAILS();

        Vec vec {{'a', 1.0}, {'b', 2.0}, {'c', 4.0}};
        double expected = 7.0;

        CHECK_EQUAL(expected, vec.NormL1());
    }

    TEST(l1_norm_calculation_with_degree_1) {
        TEST_DETAILS();

        Vec vec {{'a', 1.0}, {'b', 2.0}, {'c', 4.0}};

        // All elements in this basis have degree 1.
        double expected = 7.0;

        CHECK_EQUAL(expected, vec.NormL1(1));
    }

    TEST(l1_norm_calculation_with_degree_2) {
        TEST_DETAILS();

        Vec vec {{'a', 1.0}, {'b', 2.0}, {'c', 4.0}};

        // All elements in this basis have degree 1.
        double expected = 0.0;

        CHECK_EQUAL(expected, vec.NormL1(2));
    }

    TEST(test_add_scal_prod_key) {
        Vec v {'a', 1.0};

        v.add_scal_prod('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', 2.0}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_prod_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_prod(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', 2.0}, {'c', 2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_key) {
        Vec v {'a', 1.0};

        v.sub_scal_prod('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.sub_scal_prod(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', -2.0}, {'c', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_key) {
        Vec v {'a', 1.0};

        v.add_scal_div('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', 0.5}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_div(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', 0.5}, {'c', 0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_key) {
        Vec v {'a', 1.0};

        v.sub_scal_div('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', -0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.sub_scal_div(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', -0.5}, {'c', -0.5}};

        CHECK_EQUAL(expected, v);
    }


}


