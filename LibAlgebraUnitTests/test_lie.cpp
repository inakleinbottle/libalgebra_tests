
#include <UnitTest++/UnitTest++.h>
#include <libalgebra/libalgebra.h>
#include "time_and_details.h"

const unsigned test_alphabet_size = 5;
const unsigned test_max_degree = 3;
typedef double TestScalarType;
typedef double TestRationalType;

typedef alg::lie<TestScalarType, TestRationalType, test_alphabet_size, test_max_degree> LIE;
typedef alg::lie_basis<TestScalarType, TestRationalType, test_alphabet_size, test_max_degree> LIEBASIS;

typedef size_t LET;

class TestBasis : public alg::hall_basis<test_alphabet_size>
{
public:
    typedef double RATIONAL;
    typedef double SCALAR;


    using KEY = typename alg::hall_basis<test_alphabet_size>::KEY;
    const std::vector<std::pair<size_t, size_t>>& get_hall_set_degree_ranges() const
    {
        return hall_set_degree_ranges;
    }

    const std::vector<LET>& get_letters() const
    {
        return letters;
    }

    const std::map<LET, KEY>& get_ltk() const
    {
        return ltk;
    }

};

typedef typename TestBasis::KEY KEY;

SUITE(hall_basis_tests) {


    TEST(test_basis_created_correct_num_letters) {
        TEST_DETAILS();
        TestBasis basis {};

        CHECK_EQUAL(test_alphabet_size, basis.get_letters().size());
    }

    TEST(test_basis_created_ltk_size) {
        TEST_DETAILS();
        TestBasis basis {};

        CHECK_EQUAL(test_alphabet_size, basis.get_ltk().size());
    }


    TEST(test_growup_functionality) {
        TEST_DETAILS();
        TestBasis basis {};

        CHECK_EQUAL(2, basis.get_hall_set_degree_ranges().size());

        basis.growup(2);

        CHECK_EQUAL(3, basis.get_hall_set_degree_ranges().size());

    }

    TEST(test_degree) {
        TEST_DETAILS();
        TestBasis basis {};
        
        CHECK_EQUAL(0, basis.degree(0));
        for (KEY i=1; i<=test_alphabet_size+1; ++i) {
            CHECK_EQUAL(1, basis.degree(i));
        }
    }


    TEST(test_key_of_letter) {
        TEST_DETAILS();
        TestBasis basis {};

        for (LET i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(KEY(i), basis.keyofletter(i));
        }

    }

    TEST(test_lparent_of_key_base_letters) {
        TEST_DETAILS();
        TestBasis basis {};

        for (KEY i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(0, basis.lparent(i));
        }
        
    }

    TEST(test_rparent_of_key_base_letters) {
        TEST_DETAILS();
        TestBasis basis {};

        for (KEY i=1; i<=test_alphabet_size; ++i) {
            CHECK_EQUAL(i, basis.rparent(i));
        }
        
    }

    TEST(test_lparent_of_key_higher_key) {
        TEST_DETAILS();
        TestBasis basis {};

        KEY key = test_alphabet_size + 1;
        CHECK_EQUAL(1, basis.lparent(key));
        
    }

    TEST(test_rparent_of_key_higher_key) {
        TEST_DETAILS();
        TestBasis basis {};

        KEY key = test_alphabet_size + 1;
        CHECK_EQUAL(2, basis.rparent(key));
    }    


    TEST(test_is_letter_for_letters) {
        TEST_DETAILS();
        TestBasis basis {};

        for (KEY i=1; i<=test_alphabet_size; ++i) {
            CHECK(basis.letter(i));
        }
    }

    TEST(test_is_letter_non_letter) {
        TEST_DETAILS();
        TestBasis basis {};

        KEY key = test_alphabet_size + 1;
        CHECK(!basis.letter(key));
    }

}


class TestLieBasis : public LIEBASIS
{
public:

};

SUITE(lie_basis_tests) {

    TEST(test_product) {
        TEST_DETAILS();
        TestLieBasis basis {};

        KEY k1=1, k2=2;
        LIE result = basis.prod(k1, k2);
        
        CHECK_EQUAL(1.0, result[test_alphabet_size+1]);

    }


}






SUITE(Hall_set_size_tests) {


    struct Fixture {
        template <alg::DEG N_LETTERS, alg::DEG MAX>
        static void deg_size_test()
        {
            alg::hall_basis<N_LETTERS> hb {};
            alg::lie_basis<double, double, N_LETTERS, MAX> lb {};
                
            CHECK_EQUAL(0, lb.start_of_degree(0));
            CHECK_EQUAL(0, lb.start_of_degree(1));

            for (alg::DEG d = 2; d < MAX; ++d) {
                hb.growup(d);
                CHECK_EQUAL(hb.size(), lb.start_of_degree(d+1));
            }

        }
    };

    

    TEST_FIXTURE(Fixture, Hall_set_size_2_letters) {
        TEST_DETAILS();
        deg_size_test<2, 20>();
    }

    TEST_FIXTURE(Fixture, Hall_set_size_3_letters) {
        TEST_DETAILS();
        deg_size_test<3, 15>();
    }
    
    TEST_FIXTURE(Fixture, Hall_set_size_4_letters) {
        TEST_DETAILS();
        deg_size_test<4, 10>();
    }

    TEST_FIXTURE(Fixture, Hall_set_size_5_letters) {
        TEST_DETAILS();
        deg_size_test<5, 8>();
    }

    TEST_FIXTURE(Fixture, Hall_set_size_6_letters) {
        TEST_DETAILS();
        deg_size_test<6, 6>();
    }

    TEST_FIXTURE(Fixture, Hall_set_size_7_letters) {
        TEST_DETAILS();
        deg_size_test<7, 5>();
    }

    TEST_FIXTURE(Fixture, Hall_set_size_8_letters) {
        TEST_DETAILS();
        deg_size_test<8, 5>();
    }




}