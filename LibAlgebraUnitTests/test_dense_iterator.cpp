


#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>

#include "time_and_details.h"
#include "helpers.h"

#include <string>
#include <map>
#include <sstream>
#include <iostream>


class DITBasis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef size_t KEY;


    typedef std::map<KEY, double> MAP;

    // Default constructor
    DITBasis() {}

    KEY begin() const 
    {
        return 0;
    }

    KEY nextkey(const KEY& k) const
    {
        return k+1;
    }

    KEY end() const
    {
        return 100;
    }

    unsigned degree(const KEY &k) {
        return (k / 10) + 1;
    }

    inline static alg::DEG index_of_key(const KEY& k)
    {
        return alg::DEG{k};
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<DITBasis*, KEY> &t
    ) {
        return os << t.second;
    }

};

typedef alg::vectors::dense_vector<DITBasis, 30> Vec;



SUITE(dense_iterators) {

    TEST(test_iter) {
        TEST_DETAILS();

        Vec vec {{1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0}, {5, 5.0}};

        char key = 1;
        double coeff = 1.0;
        for (auto it : vec) {
            CHECK_EQUAL(key, it.first);
            CHECK_EQUAL(coeff, it.second);
            ++key;
            coeff += 1.0;
        }

    }

    TEST(test_iterator_with_zeros_empty) {
        TEST_DETAILS();
        Vec vec;

        CHECK(vec.begin() == vec.end());
    }

    TEST(test_iterator_one_value) {
        TEST_DETAILS();
        Vec vec {15, 1.0};
        typename Vec::iterator itbegin, itend;
        itbegin = vec.begin();
        itend = vec.end();

        REQUIRE(itbegin != itend);
        CHECK_EQUAL(15, itbegin->first);
        CHECK_EQUAL(1.0, itbegin->second);

        ++itbegin;
        CHECK(itend == itbegin);
    }

    TEST(test_iterator_two_values) {
        TEST_DETAILS();
        Vec vec {{10, 1.0}, {20, 2.0}};
        typename Vec::iterator itbegin, itend;
        itbegin = vec.begin();
        itend = vec.end();

        REQUIRE(itbegin != itend);
        CHECK_EQUAL(10, itbegin->first);
        CHECK_EQUAL(1.0, itbegin->second);

        ++itbegin;
        REQUIRE(itbegin != itend);
        CHECK_EQUAL(20, itbegin->first);
        CHECK_EQUAL(2.0, itbegin->second);

        ++itbegin;
        CHECK(itend == itbegin);

    }

    TEST(test_iterate_over_specific_degree) {
        TEST_DETAILS();
        Vec vec {
            {1, 1.0}, {2, 2.0}, {3, 3.0},  // deg 1
            {11, 1.1}, {12, 2.1}, {13, 3.1},  // deg 2
            {21, 1.2}, {22, 2.2}, {23, 3.2}   // deg 3
        };

        size_t degree = 2;
        typename Vec::iterator itbegin, itend;
        itbegin = vec.deg_begin(degree);
        itend = vec.deg_end(degree);

        REQUIRE(itbegin != vec.begin());
        REQUIRE(itend != vec.end());
        REQUIRE(itbegin != itend);

        size_t i=0;
        for(; itbegin != itend; ++itbegin) {
            ++i;
            CHECK_EQUAL(10+i, itbegin->first);
            CHECK_EQUAL(i+0.1, itbegin->second);
        }
    }


    TEST(test_dense_iterator_modify_value) {
        TEST_DETAILS();

        Vec vec {
            {1, 1.0}, {2, 2.0}, {3, 3.0},  // deg 1
            {11, 1.1}, {12, 2.1}, {13, 3.1},  // deg 2
            {21, 1.2}, {22, 2.2}, {23, 3.2}   // deg 3
        };

        Vec expected = {
            {2, 1.0}, {3, 2.0},  // deg 1
            {11, 0.1}, {12, 1.1}, {13, 2.1},  // deg 2
            {21, 0.2}, {22, 1.2}, {23, 2.2}   // deg 3
        };

        for (typename Vec::iterator it=vec.begin(); it!=vec.end(); ++it)
            it->second -= 1;
        
        CHECK_VEC_CLOSE(expected, vec, 2.0e-15);


    }

}