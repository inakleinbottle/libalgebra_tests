
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
    typedef double SCALAR;
    typedef std::string KEY;
    typedef double SCA;
    typedef std::map<KEY, SCA> MAP;

    typedef alg::algebra<Basis, 26> ALG;
    typedef std::pair<KEY, KEY> CACHE_KEY;
    typedef std::map<CACHE_KEY, ALG> CACHE; 

    // Default constructor
    Basis() : _cache{} {
    }

    unsigned degree(const KEY &k) {
        return k.size();
    }

    inline KEY begin()
    {
        return "";
    }

    inline KEY end()
    {
        return "zzz";
    }

    inline KEY nextkey(const KEY &key)
    {
        if (key.empty())
            return "a";

        return std::string { char(key.back() + 1) };
        
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


    inline static bool comp(const KEY& k1, const KEY& k2)
    {
        if (k1.size() < k2.size())
            return true;
        else if (k1.size() > k2.size())
            return false;
        
        for (int i=0; i<k1.size(); ++i) {
            if (k1[i] < k2[i])
                return true;
            else if (k1[i] > k2[i])
                return false;
        }
        return true;
    }

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
        Alg actual = a*b;

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

    TEST(test_inplace_multiplication_reverse_order) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg b {"b", 2.0};
        
        Alg expected {"ba", 2.0};
        b *= a;

        CHECK_EQUAL(expected, b);

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

    TEST(test_truncate_by_degree_upper_bound) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        a["b"] = 1.0;
        a["ab"] = 1.0;
        a["aba"] = 1.0;

        Alg expected {"a", 1.0};
        expected["b"] = 1.0;

        Alg c = a.truncate(0, 1);

        CHECK_EQUAL(expected, c);
    }

    TEST(test_truncate_by_degree_lower_bound) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        a["b"] = 1.0;
        a["ab"] = 1.0;
        a["aba"] = 1.0;

        Alg expected {"ab", 1.0};
        expected["aba"] = 1.0;

        Alg c = a.truncate(2, 5);

        CHECK_EQUAL(expected, c);
    }


    TEST(test_truncate_by_degree_both_bounds) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        a["b"] = 1.0;
        a["ab"] = 1.0;
        a["aba"] = 1.0;

        Alg expected {"ab", 1.0};

        Alg c = a.truncate(2, 2);

        CHECK_EQUAL(expected, c);
    }

    TEST(test_commutator_formulation) {
        TEST_DETAILS();
        Alg a {"a", 1.0};
        Alg b {"b", 1.0};

        Alg c = commutator(a, b);

        Alg expected {"ab", 1.0};
        expected["ba"] = -1.0;

        CHECK_EQUAL(expected, c);

    }

    TEST(test_add_mul) {
        TEST_DETAILS();

        Alg a {"a", 1.0};
        Alg b {"b", 1.0};
        Alg c {"c", 2.0};

        Alg result = a.add_mul(b, c);

        Alg expected {"a", 1.0};
        expected["bc"] = 2.0;

        CHECK_EQUAL(expected, result);
    }

    TEST(test_sub_mul) {
        TEST_DETAILS();

        Alg a {"a", 1.0};
        Alg b {"b", 1.0};
        Alg c {"c", 2.0};

        Alg result = a.sub_mul(b, c);

        Alg expected {"a", 1.0};
        expected["bc"] = -2.0;

        CHECK_EQUAL(expected, result);
    }

/*    
    //This test causes a compile error on line 333.
   //There seems to be a missing (unused) argument of type Identity<0>

    TEST(test_mul_scal_prod) {
        TEST_DETAILS();

        Alg a {"a", 1.0};
        Alg b {"b", 1.0};
        double sca = 2.0;

        Alg result = a.mul_scal_prod(b, sca);

        Alg expected {"ab", 2.0};

        CHECK_EQUAL(expected, result);
    }
    */

    TEST(test_mul_scal_div) {
        TEST_DETAILS();

        Alg a {"a", 1.0};
        Alg b {"b", 2.0};
        double sca = 2.0;

        Alg result = a.mul_scal_div(b, sca);

        Alg expected {"ab", 1.0};

        CHECK_EQUAL(expected, result);
    }

} 