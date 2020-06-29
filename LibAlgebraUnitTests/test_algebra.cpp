
#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>
#include "time_and_details.h"

#include <map>
#include <cstddef>
#include <string>


// Minimal algebra basis class
class Basis {
public:
    typedef double RATIONAL;
    typedef std::string KEY;
    typedef double SCA;
    typedef std::map<KEY, SCA> MAP;



    // Default constructor
    Basis() : _cache{} {}

    unsigned degree(const KEY &k) {
        return k.size();
    }


    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<Basis*, KEY> &t
    ) {
        return os << t.second;
    }

    // All the above are necessary for the sparse_vector basis.

    // Max degree is ignored for the purposes of testing here,
    // but is necessary for instantiating an algebra object.
    static const unsigned MAX_DEGREE = 3; 

    typedef alg::algebra<Basis> ALG;
    typedef std::pair<KEY, KEY> CACHE_KEY;
    typedef std::map<CACHE_KEY, ALG> CACHE; 

    // This is a very basic product function that uses a cache
    // to store the results of products of two basis vectors.
    inline const ALG& prod(const KEY& k1, const KEY& k2) {
        typename CACHE::iterator it;

        CACHE_KEY ck (k1, k2);
        it = _cache.find(ck);
        if (it == _cache.end()) {
            _cache[ck] = ALG{k1 + k2};
            return _cache[ck];            
        } else {
            return it->second;
        }

    }
    
private:
    CACHE _cache;

};

typedef typename Basis::ALG Alg;


SUITE(test_algebra){ 

    TEST(test_binary_multiplication_different_elements) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg b {"b", 2.0};

        Alg expected {"ab", 2.0};

        CHECK_EQUAL(expected, a * b);
    }

    TEST(test_binary_multiplication_zero_element) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg e {};

        CHECK_EQUAL(e, e*a);
        CHECK_EQUAL(e, a*e);
    }

    TEST(test_inplace_multiplication_different_elements) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg b {"b", 2.0};

        Alg expected {"ab", 2.0};
        a *= b;

        CHECK_EQUAL(expected, a);
    }

    TEST(test_inplace_left_multiplication_zero_element) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg e {};

        e *= a;
        Alg expected {};

        CHECK_EQUAL(expected, e);
    }

    TEST(test_inplace_right_multiplication_zero_element) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg e {};

        e *= a;
        Alg expected {};

        CHECK_EQUAL(expected, e);
    }


    TEST(test_degree_single_item) {
        TEST_DETAILS();
        Alg a {"a", 1.0};

        CHECK_EQUAL(1, a.degree());
    }

    TEST(test_degree_multiple_items) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        a["b"] = 1.0;
        a["ab"] = 1.0;

        CHECK_EQUAL(2, a.degree());
    }



}