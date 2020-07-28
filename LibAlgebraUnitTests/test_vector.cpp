
#include <UnitTest++/UnitTest++.h>
#include <libalgebra/vectors/vectors.h>


#include "time_and_details.h"



#include <map>
#include <utility>
#include <iostream>





class VTBasis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef char KEY;
    typedef std::map<char, double> MAP;

    // Default constructor
    VTBasis() {}

    KEY begin() const 
    {
        return 'a';
    }

    KEY nextkey(const KEY& k) const
    {
        return char(k + 1);
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
        const std::pair<VTBasis*, KEY> &t
    ) {
        return os << t.second;
    }

    inline static bool comp(const KEY& k1, const KEY& k2)
	{
		return (k1 == k2) || k1 < k2;
	}

};

typedef alg::vectors::vector<VTBasis, 2> Vec;


SUITE(vector_tests) {

    TEST(test_vector_setup_keys) {
        TEST_DETAILS();
        Vec v {};

        const std::vector<char>& keys = v.get_dense_keys();
        
        REQUIRE(2, keys.size());
        CHECK_EQUAL('a', keys[0]);
        CHECK_EQUAL('b', keys[1]);
    }

    TEST(test_vector_setup_coeffs) {
        TEST_DETAILS();
        Vec v {};

        const std::vector<double>& coeffs = v.get_dense_coeffs();

        CHECK_EQUAL(0.0, coeffs[0]);
        CHECK_EQUAL(0.0, coeffs[1]);
    }

    TEST(test_initilizer_list_creation) {
        TEST_DETAILS();
        Vec v {{'a', 1.0}, {'b', 2.0}, {'c', 3.0}};

        CHECK_EQUAL(1.0, v['a']);
        CHECK_EQUAL(2.0, v['b']);
        CHECK_EQUAL(3.0, v['c']);
    }

    TEST(test_element_access_dense_part) {
        TEST_DETAILS();
        Vec v {'a', 1.0};

        CHECK_EQUAL(1.0, v['a']);
        CHECK_EQUAL(0.0, v['b']);
    }

    TEST(test_element_access_sparse_part) {
        TEST_DETAILS();
        Vec v {'c', 1.0};

        CHECK_EQUAL(1.0, v['c']);
    }

    TEST(test_inserting_coordinate_dense_part) {
        TEST_DETAILS();
        Vec vec {};
        
        vec['a'] = 1.0;
        try {
            CHECK(vec.dense_part()['a'] == 1.0);
        }catch (alg::vectors::KeyNotFoundError) {
            std::cerr << "Uncaught key not found error" << "\n";
        }
    }

    TEST(test_inserting_coordinate_sparse_part) {
        TEST_DETAILS();
        Vec vec {};
        
        vec['c'] = 1.0;
        CHECK(vec.sparse_part()['c'] == 1.0);   
    }

   TEST(test_inplace_scalar_multiplication) {
        TEST_DETAILS();
        Vec vec {{'a', 1.0}, {'c', 1.0}}; 
        double scalar = 2.0;

        vec *= scalar;
        CHECK_EQUAL(2.0, vec['a']);
        CHECK_EQUAL(2.0, vec['c']);
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
        Vec vec {{'a', 1.0}, {'b', 1.0}};
        double scalar = 2.0;

        Vec nvec = vec*scalar;
        Vec expected {{'a', 2.0}, {'b', 2.0}};

        CHECK_EQUAL(expected, nvec);
    }
       
    TEST(test_derived_scalar_multiplication_neutral_element) {
        TEST_DETAILS();
        Vec neut {};
        double scalar = 2.0;

        Vec nvec = neut * scalar;

        CHECK_EQUAL(neut, nvec);
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
        Vec vec1 {{'a', 1.0}, {'c', 1.0}};
        Vec vec2 {{'a', 2.0}, {'c', 1.0}};

        Vec nvec = vec1 + vec2;
        Vec expected {{'a', 3.0}, {'c', 2.0}};

        CHECK_EQUAL(expected, nvec);
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

    TEST(test_equality_equal_vectors_dense_part) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        CHECK_EQUAL(vec1, vec2);
    }

    TEST(test_equality_not_equal_keys_dense_part) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'b', 1.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_not_equal_values_dense_part) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 2.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_equal_vectors_sparse_part) {
        TEST_DETAILS();
        Vec vec1 {'a', 1.0};
        Vec vec2 {'a', 1.0};

        CHECK_EQUAL(vec1, vec2);
    }

    TEST(test_equality_not_equal_keys_sparse_part) {
        TEST_DETAILS();
        Vec vec1 {'c', 1.0};
        Vec vec2 {'d', 1.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_not_equal_values_sparse_part) {
        TEST_DETAILS();
        Vec vec1 {'c', 1.0};
        Vec vec2 {'c', 2.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_equal_sparse_different_dense) {
        TEST_DETAILS();
        Vec v1 {{'a', 1.0}, {'c', 1.0}};
        Vec v2 {{'b', 2.0}, {'c', 1.0}};

        CHECK(!(v1 == v2));
    }

    TEST(test_equality_equal_dense_different_sparse) {
        TEST_DETAILS();
        Vec v1 {{'a', 1.0}, {'c', 1.0}};
        Vec v2 {{'a', 1.0}, {'d', 2.0}};

        CHECK(!(v1 == v2));
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


    TEST(test_iterator) {
        TEST_DETAILS();

        Vec vec {{'a', 1.0}, {'b', 2.0}, {'c', 3.0}, {'d', 4.0}};

        typename Vec::iterator it = vec.begin();
        double c = 1.0;
        char key = 'a';
        for (; it!=vec.end(); ++it) {
            CHECK_EQUAL(key, it->first);
            CHECK_EQUAL(c, it->second);
            c += 1.0;
            ++key;
        }

    }

    TEST(test_add_scal_prod_dense_key) {
        Vec v {'a', 1.0};

        v.add_scal_prod('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', 2.0}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_prod_sparse_key) {
        Vec v {'a', 1.0};

        v.add_scal_prod('c', 2.0);
        Vec expected {{'a', 1.0}, {'c', 2.0}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_prod_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_prod(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', 2.0}, {'c', 2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_dense_key) {
        Vec v {'a', 1.0};

        v.sub_scal_prod('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_sparse_key) {
        Vec v {'a', 1.0};

        v.sub_scal_prod('c', 2.0);
        Vec expected {{'a', 1.0}, {'c', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.sub_scal_prod(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', -2.0}, {'c', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_dense_key) {
        Vec v {'a', 1.0};

        v.add_scal_div('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', 0.5}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_sparse_key) {
        Vec v {'a', 1.0};

        v.add_scal_div('c', 2.0);
        Vec expected {{'a', 1.0}, {'c', 0.5}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_vector) {
        Vec v {'a', 1.0};
        Vec rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_div(rhs, 2.0);
        Vec expected {{'a', 1.0}, {'b', 0.5}, {'c', 0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_dense_key) {
        Vec v {'a', 1.0};

        v.sub_scal_div('b', 2.0);
        Vec expected {{'a', 1.0}, {'b', -0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_sparse_key) {
        Vec v {'a', 1.0};

        v.sub_scal_div('c', 2.0);
        Vec expected {{'a', 1.0}, {'c', -0.5}};

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




typedef alg::vectors::vector<VTBasis, 0> Vec0;


SUITE(vector_transition_0_tests) {

    TEST(test_inserting_coordinate) {
        TEST_DETAILS();
        Vec0 vec {};
        
        vec['a'] = 1.0;
        CHECK(vec['a'] == 1.0);   
    }

    TEST(text_formatting) {        
        TEST_DETAILS();
        std::basic_ostringstream<char> stream {};
        Vec0 vec {'a', 1.0};

        stream << vec;

        CHECK_EQUAL("{ 1(a) }", stream.str());
    }

    TEST(test_unary_minus) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        vec['b'] = 2.0;

        Vec0 nvec = -vec;
        CHECK_EQUAL(-1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_inplace_scalar_multiplication) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        double scalar = 2.0;

        vec *= scalar;
        CHECK_EQUAL(2.0, vec['a']);
    }    
    
    TEST(test_inplace_scalar_multiplication_neutral_element) {
        TEST_DETAILS();
        Vec0 neut {};
        double scalar = 2.0;

        neut *= scalar;

        CHECK_EQUAL(Vec0 {}, neut);
    }

    TEST(test_derived_scalar_multiplication) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        double scalar = 2.0;

        Vec0 nvec = vec*scalar;

        CHECK_EQUAL(2.0, nvec['a']);
    }

        
    TEST(test_derived_scalar_multiplication_neutral_element) {
        TEST_DETAILS();
        Vec0 neut {};
        double scalar = 2.0;

        Vec0 nvec = neut * scalar;

        CHECK_EQUAL(neut, nvec);
    }


    TEST(test_inplace_addition_same_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(3.0, vec1['a']);
    }

    TEST(test_inplace_addition_different_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'b', 2.0};

        vec1 += vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(2.0, vec1['b']);
    }

    TEST(test_inplace_addition_interlaced_order) {
        TEST_DETAILS();
        Vec0 vec1 {'b', 1.0};
        vec1['d'] = 1.0;
        Vec0 vec2 {'a', 1.0};
        vec2['c'] = 1.0;

        Vec0 expected {'a', 1.0};
        expected['b'] = 1.0;
        expected['c'] = 1.0;
        expected['d'] = 1.0;

        vec1 += vec2;
        CHECK_EQUAL(expected, vec1);
    }

    TEST(test_inplace_addition_neutral_element) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        Vec0 neut {};

        vec += neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_binary_addition_same_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 2.0};

        Vec0 nvec = vec1 + vec2;

        CHECK_EQUAL(3.0, nvec['a']);
    }

    TEST(test_derived_addition_different_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'b', 2.0};

        Vec0 nvec = vec1 + vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(2.0, nvec['b']);
    }

    TEST(test_binary_addition_interlaced_order) {
        TEST_DETAILS();
        Vec0 vec1 {'b', 1.0};
        vec1['d'] = 1.0;
        Vec0 vec2 {'a', 1.0};
        vec2['c'] = 1.0;

        Vec0 expected {'a', 1.0};
        expected['b'] = 1.0;
        expected['c'] = 1.0;
        expected['d'] = 1.0;

        Vec0 nvec = vec1 + vec2;
        CHECK_EQUAL(expected, nvec);
    }
        
    TEST(test_derived_addition_neutral_element) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        Vec0 neut {};

        Vec0 nvec = vec + neut;
        CHECK_EQUAL(1.0, nvec['a']);
    }

    TEST(test_inplace_subtraction_same_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(-1.0, vec1['a']);
    }

    TEST(test_inplace_subtract_same_keys_to_0) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 1.0};

        vec1 -= vec2;

        CHECK_EQUAL(0.0, vec1['a']);
    }

    TEST(test_inplace_subtraction_different_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'b', 2.0};

        vec1 -= vec2;

        CHECK_EQUAL(1.0, vec1['a']);
        CHECK_EQUAL(-2.0, vec1['b']);
    }

    TEST(test_inplace_subtraction_neutral_element) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        Vec0 neut {};

        vec -= neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_derived_subtraction_same_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 2.0};

        Vec0 nvec = vec1 - vec2;

        CHECK_EQUAL(-1.0, nvec['a']);
    }

    TEST(test_derived_subtraction_same_keys_to_0) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 1.0};

        Vec0 nvec = vec1 - vec2;

        CHECK_EQUAL(0.0, nvec['a']);
    }

    TEST(test_derived_subtraction_different_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'b', 2.0};

        Vec0 nvec = vec1 - vec2;

        CHECK_EQUAL(1.0, nvec['a']);
        CHECK_EQUAL(-2.0, nvec['b']);
    }

    TEST(test_derived_subtraction_neutral_element) {
        TEST_DETAILS();
        Vec0 vec {'a', 1.0};
        Vec0 neut {};

        Vec0 nvec = vec - neut;
        CHECK_EQUAL(1.0, vec['a']);
    }

    TEST(test_equality_equal_vectors) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 1.0};

        CHECK_EQUAL(vec1, vec2);
    }

    TEST(test_equality_not_equal_keys) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'b', 1.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_not_equal_values) {
        TEST_DETAILS();
        Vec0 vec1 {'a', 1.0};
        Vec0 vec2 {'a', 2.0};

        CHECK(!(vec1 == vec2));
    }

    TEST(test_equality_neutral_element) {
        TEST_DETAILS();
        Vec0 neut1 {}, neut2 {};

        CHECK_EQUAL(neut1, neut2);
    }
    
    TEST(l1_norm_calculation) {
        TEST_DETAILS();

        Vec0 vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;
        double expected = 6.0;

        CHECK_EQUAL(expected, vec.NormL1());
    }

    TEST(l1_norm_calculation_with_degree_1) {
        TEST_DETAILS();

        Vec0 vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;

        // All elements in this basis have degree 1.
        double expected = 6.0;

        CHECK_EQUAL(expected, vec.NormL1(1));
    }

    TEST(l1_norm_calculation_with_degree_2) {
        TEST_DETAILS();

        Vec0 vec {'a', 1.0};
        vec['b'] = 2.0;
        vec['c'] = 3.0;

        // All elements in this basis have degree 1.
        double expected = 0.0;

        CHECK_EQUAL(expected, vec.NormL1(2));
    }

        TEST(test_add_scal_prod_key) {
        Vec0 v {'a', 1.0};

        v.add_scal_prod('b', 2.0);
        Vec0 expected {{'a', 1.0}, {'b', 2.0}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_prod_vector) {
        Vec0 v {'a', 1.0};
        Vec0 rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_prod(rhs, 2.0);
        Vec0 expected {{'a', 1.0}, {'b', 2.0}, {'c', 2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_key) {
        Vec0 v {'a', 1.0};

        v.sub_scal_prod('b', 2.0);
        Vec0 expected {{'a', 1.0}, {'b', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_prod_vector) {
        Vec0 v {'a', 1.0};
        Vec0 rhs {{'b', 1.0}, {'c', 1.0}};

        v.sub_scal_prod(rhs, 2.0);
        Vec0 expected {{'a', 1.0}, {'b', -2.0}, {'c', -2.0}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_key) {
        Vec0 v {'a', 1.0};

        v.add_scal_div('b', 2.0);
        Vec0 expected {{'a', 1.0}, {'b', 0.5}};
        
        CHECK_EQUAL(expected, v);
    }

    TEST(test_add_scal_div_vector) {
        Vec0 v {'a', 1.0};
        Vec0 rhs {{'b', 1.0}, {'c', 1.0}};

        v.add_scal_div(rhs, 2.0);
        Vec0 expected {{'a', 1.0}, {'b', 0.5}, {'c', 0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_key) {
        Vec0 v {'a', 1.0};

        v.sub_scal_div('b', 2.0);
        Vec0 expected {{'a', 1.0}, {'b', -0.5}};

        CHECK_EQUAL(expected, v);
    }

    TEST(test_sub_scal_div_vector) {
        Vec0 v {'a', 1.0};
        Vec0 rhs {{'b', 1.0}, {'c', 1.0}};

        v.sub_scal_div(rhs, 2.0);
        Vec0 expected {{'a', 1.0}, {'b', -0.5}, {'c', -0.5}};

        CHECK_EQUAL(expected, v);
    }


}



















