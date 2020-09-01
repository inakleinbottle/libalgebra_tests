#pragma once


#include <addons/gmpwrapper.h>


namespace {

    template<bool>
    struct check_vec_almost_equal_impl {
        template<class VECT>
        inline static void check(const VECT& expected, const VECT& actual, double tol)
        {
            CHECK_EQUAL(expected, actual);
        }
    };

    template<>
    struct check_vec_almost_equal_impl<false> {
        template<class VECT>
        inline static void check(const VECT& expected, const VECT& actual, double tol)
        {
            VECT err = actual - expected;
            for (auto& elt : err) {
                CHECK_CLOSE(0.0, elt.second, tol);
            }
        }
    };

    template <class VECT>
    inline void check_vec_almost_equal(const VECT& expected, const VECT& actual, double tol)
    {
        typedef typename VECT::SCALAR S;
        check_vec_almost_equal_impl<std::is_same<S, mpq_class>::value>::
            check(expected, actual, tol);
    }


}

#define CHECK_VEC_CLOSE(expected, actual, tol) check_vec_almost_equal(expected, actual, tol)