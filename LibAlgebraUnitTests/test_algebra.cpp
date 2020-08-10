
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

    typedef alg::algebra<Basis, 12> ALG;
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
        return "z";
    }

    inline KEY nextkey(const KEY &key)
    {
        if (key.empty())
            return "a";

        return std::string { char(key.back() + 1) };
        
    }

    inline static alg::DEG index_of_key(const KEY& k)
    {
        if (k == "") return 0;
        return alg::DEG{k[0] - 'a' + 1};
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

    TEST(test_initializer_list_creation) {
        Alg a = {
            {"a", 1.0}, {"b", 1.0}, {"c", 1.0},
            {"m", 1.0}, {"n", 1.0}, {"o", 1.0}
        };

        CHECK_EQUAL(6, a.size());
        
        CHECK_EQUAL(1.0, a["a"]);
        CHECK_EQUAL(1.0, a["b"]);
        CHECK_EQUAL(1.0, a["c"]);
        CHECK_EQUAL(1.0, a["m"]);
        CHECK_EQUAL(1.0, a["n"]);
        CHECK_EQUAL(1.0, a["o"]);

    }

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





SUITE(algebra_large) {

    struct VectorFixture {
        Alg alg1, alg2;

        VectorFixture()
        {
            alg1 = {
                {"b", 0.3}, {"c", -0.1}, {"d", -0.6},
                {"f", 1.1}, {"i", -0.7}, {"j", 1.7},
                {"k", 1.7}, {"m", 0.4}, {"o", 0.3},
                {"p", -0.7}, {"r", -0.0}, {"x", -0.0},
                
            };
            alg2 = {
                {"c", 0.5}, {"e", 0.5}, {"f", -0.1},
                {"g", 0.1}, {"h", 0.2}, {"i", 0.6},
                {"j", -0.3}, {"l", -0.4}, {"o", 0.3},
                {"p", 1.9}, {"r", -0.2}, {"t", -0.2},
                {"u", -0.7}, {"w", 0.7}, {"y", 1.1},
                {"z", 1.6}, 
            };
        }
    };

    TEST_FIXTURE(VectorFixture, test_inplace_mul) {
        TEST_DETAILS();
        Alg expected = {
        {"bc", 0.15}, {"be", 0.15}, {"bf", -0.03},
            {"bg", 0.03}, {"bh", 0.06}, {"bi", 0.18},
            {"bj", -0.09}, {"bl", -0.12}, {"bo", 0.09},
            {"bp", 0.57}, {"br", -0.06}, {"bt", -0.06},
            {"bu", -0.21}, {"bw", 0.21}, {"by", 0.33},
            {"bz", 0.48}, {"cc", -0.05}, {"ce", -0.05},
            {"cf", 0.010000000000000002}, {"cg", -0.010000000000000002}, {"ch", -0.020000000000000004},
            {"ci", -0.06}, {"cj", 0.03}, {"cl", 0.04000000000000001},
            {"co", -0.03}, {"cp", -0.19}, {"cr", 0.020000000000000004},
            {"ct", 0.020000000000000004}, {"cu", 0.06999999999999999}, {"cw", -0.06999999999999999},
            {"cy", -0.11000000000000001}, {"cz", -0.16000000000000003}, {"dc", -0.3},
            {"de", -0.3}, {"df", 0.06}, {"dg", -0.06},
            {"dh", -0.12}, {"di", -0.36}, {"dj", 0.18},
            {"dl", 0.24}, {"do", -0.18}, {"dp", -1.14},
            {"dr", 0.12}, {"dt", 0.12}, {"du", 0.42},
            {"dw", -0.42}, {"dy", -0.66}, {"dz", -0.96},
            {"fc", 0.55}, {"fe", 0.55}, {"ff", -0.11000000000000001},
            {"fg", 0.11000000000000001}, {"fh", 0.22000000000000003}, {"fi", 0.66},
            {"fj", -0.33}, {"fl", -0.44000000000000006}, {"fo", 0.33},
            {"fp", 2.09}, {"fr", -0.22000000000000003}, {"ft", -0.22000000000000003},
            {"fu", -0.77}, {"fw", 0.77}, {"fy", 1.2100000000000002},
            {"fz", 1.7600000000000002}, {"ic", -0.35}, {"ie", -0.35},
            {"if", 0.06999999999999999}, {"ig", -0.06999999999999999}, {"ih", -0.13999999999999999},
            {"ii", -0.42}, {"ij", 0.21}, {"il", 0.27999999999999997},
            {"io", -0.21}, {"ip", -1.3299999999999998}, {"ir", 0.13999999999999999},
            {"it", 0.13999999999999999}, {"iu", 0.48999999999999994}, {"iw", -0.48999999999999994},
            {"iy", -0.77}, {"iz", -1.1199999999999999}, {"jc", 0.85},
            {"je", 0.85}, {"jf", -0.17}, {"jg", 0.17},
            {"jh", 0.34}, {"ji", 1.02}, {"jj", -0.51},
            {"jl", -0.68}, {"jo", 0.51}, {"jp", 3.23},
            {"jr", -0.34}, {"jt", -0.34}, {"ju", -1.19},
            {"jw", 1.19}, {"jy", 1.87}, {"jz", 2.72},
            {"kc", 0.85}, {"ke", 0.85}, {"kf", -0.17},
            {"kg", 0.17}, {"kh", 0.34}, {"ki", 1.02},
            {"kj", -0.51}, {"kl", -0.68}, {"ko", 0.51},
            {"kp", 3.23}, {"kr", -0.34}, {"kt", -0.34},
            {"ku", -1.19}, {"kw", 1.19}, {"ky", 1.87},
            {"kz", 2.72}, {"mc", 0.2}, {"me", 0.2},
            {"mf", -0.04000000000000001}, {"mg", 0.04000000000000001}, {"mh", 0.08000000000000002},
            {"mi", 0.24}, {"mj", -0.12}, {"ml", -0.16000000000000003},
            {"mo", 0.12}, {"mp", 0.76}, {"mr", -0.08000000000000002},
            {"mt", -0.08000000000000002}, {"mu", -0.27999999999999997}, {"mw", 0.27999999999999997},
            {"my", 0.44000000000000006}, {"mz", 0.6400000000000001}, {"oc", 0.15},
            {"oe", 0.15}, {"of", -0.03}, {"og", 0.03},
            {"oh", 0.06}, {"oi", 0.18}, {"oj", -0.09},
            {"ol", -0.12}, {"oo", 0.09}, {"op", 0.57},
            {"or", -0.06}, {"ot", -0.06}, {"ou", -0.21},
            {"ow", 0.21}, {"oy", 0.33}, {"oz", 0.48},
            {"pc", -0.35}, {"pe", -0.35}, {"pf", 0.06999999999999999},
            {"pg", -0.06999999999999999}, {"ph", -0.13999999999999999}, {"pi", -0.42},
            {"pj", 0.21}, {"pl", 0.27999999999999997}, {"po", -0.21},
            {"pp", -1.3299999999999998}, {"pr", 0.13999999999999999}, {"pt", 0.13999999999999999},
            {"pu", 0.48999999999999994}, {"pw", -0.48999999999999994}, {"py", -0.77},
            {"pz", -1.1199999999999999}, {"rc", 0.0}, {"re", 0.0},
            {"rf", 0.0}, {"rg", 0.0}, {"rh", 0.0},
            {"ri", 0.0}, {"rj", 0.0}, {"rl", 0.0},
            {"ro", 0.0}, {"rp", 0.0}, {"rr", 0.0},
            {"rt", 0.0}, {"ru", 0.0}, {"rw", 0.0},
            {"ry", 0.0}, {"rz", 0.0}, {"xc", 0.0},
            {"xe", 0.0}, {"xf", 0.0}, {"xg", 0.0},
            {"xh", 0.0}, {"xi", 0.0}, {"xj", 0.0},
            {"xl", 0.0}, {"xo", 0.0}, {"xp", 0.0},
            {"xr", 0.0}, {"xt", 0.0}, {"xu", 0.0},
            {"xw", 0.0}, {"xy", 0.0}, {"xz", 0.0},
        };

        alg1 *= alg2;

        CHECK_EQUAL(expected, alg1);
    }

}