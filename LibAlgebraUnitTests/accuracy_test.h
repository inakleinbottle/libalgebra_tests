//
// Created by sam on 22/01/2021.
//

#ifndef LIBALGEBRAUNITTESTS_ACCURACY_TEST_H
#define LIBALGEBRAUNITTESTS_ACCURACY_TEST_H

#include <random>

#include <UnitTest++/UnitTest++.h>


#include "alg_framework.h"


template <alg::DEG Width, alg::DEG Depth>
class AccuracyTest {
    using RationalF = alg_framework<Depth, Width, Rational>;
    using FloatF = alg_framework<Depth, Width, SPReal>;
    using DoubleF = alg_framework<Depth, Width, DPReal>;


public:

    typedef typename FloatF::DEG DEG;
    typedef typename alg::DIMN DIMN;
    typedef typename FloatF::LET LET;

    typedef typename RationalF::S Sca;
    typedef typename RationalF::RAT Rat;

    typedef typename RationalF::TENSOR::KEY TKey;


    struct uniform
    {
        uniform() : upper{1}, lower{0}, diff{1}
        {}

        uniform(Sca lower, Sca upper) : upper{upper}, lower{lower}, diff{upper - lower}
        {
            assert(upper > lower);
        }

        template <typename RNG>
        inline Sca rand(RNG rng)
        {
            return Sca(rng()) / Sca(rng.max());
        }

        template <typename RNG>
        inline Sca operator()(RNG &rng)
        {
            return rand(rng)*diff + lower;
        }

    private:
        Sca lower, upper;
        Sca diff;
    };



#define DECLARE_BINARY_OPERATOR(RT, OP, OT)                               \
    inline RT operator OP(const OT& other) const                          \
    {                                                                     \
        return RT(                                                        \
            m_rdata OP other.m_rdata,                                     \
            m_fdata OP other.m_fdata,                                     \
            m_ddata OP other.m_ddata                                      \
            );                                                            \
    }

#define DECLARE_INPLACE_BINARY_OP(RT, OP, OT)                             \
    inline RT& operator OP(const OT& other)                               \
    {                                                                     \
        m_rdata OP other.m_rdata;                                         \
        m_fdata OP other.m_fdata;                                         \
        m_ddata OP other.m_ddata;                                         \
        return *this;                                                     \
    }

#define DECLARE_UMINUS(RT)                                                \
    inline RT operator-() const                                           \
    {                                                                     \
        return RT(-m_rdata, -m_fdata. -m_ddata);                          \
    }

    class Maps;

    class Lie
    {
        typedef typename RationalF::LIE RLie;
        typedef typename FloatF::LIE FLie;
        typedef typename DoubleF::LIE DLie;

        RLie m_rdata;
        FLie m_fdata;
        DLie m_ddata;

        Lie(RLie r, FLie f, DLie d) : m_rdata(r), m_fdata(f), m_ddata(d)
        {}

        friend class Maps;

    public:

        Lie() : m_rdata(), m_fdata(), m_ddata()
        {}

        template<typename F>
        void fill_with(F fn, DEG depth)
        {
            using VBP = alg::vectors::VectorBasisProperties<typename RLie::BASIS>;
            DIMN size = VBP::DegProp::start_of_degree(depth+1);
            std::vector<Sca> tmp;
            tmp.reserve(size);


            auto basis = RLie::basis;
            typename RLie::KEY key {basis.begin()};
            for (DIMN i=0; i<size; ++i) {
                Sca val = fn(i);
                m_rdata[key] = val;
                key = basis.nextkey(key);
                tmp.push_back(val);
            }

            REQUIRE CHECK(tmp.size() == size);

            m_fdata.fill_with([&] (DIMN i) { return float(tmp[i]); }, depth);
            m_ddata.fill_with([&] (DIMN i) { return double(tmp[i]); }, depth);
        }


//        DECLARE_UMINUS(Tensor)

        DECLARE_BINARY_OPERATOR(Lie, +, Lie)
        DECLARE_BINARY_OPERATOR(Lie, -, Lie)

        DECLARE_BINARY_OPERATOR(Lie, *, Sca)
        DECLARE_BINARY_OPERATOR(Lie, /, Rat)

        DECLARE_BINARY_OPERATOR(Lie, *, Lie)

        DECLARE_INPLACE_BINARY_OP(Lie, +=, Lie)
        DECLARE_INPLACE_BINARY_OP(Lie, -=, Lie)

        DECLARE_INPLACE_BINARY_OP(Lie, +=, Sca)
        DECLARE_INPLACE_BINARY_OP(Lie, /=, Rat)

        DECLARE_INPLACE_BINARY_OP(Lie, *=, Lie)

    };

    class Tensor
    {
        typedef typename RationalF::TENSOR RTensor;
        typedef typename FloatF::TENSOR FTensor;
        typedef typename DoubleF::TENSOR DTensor;

        RTensor m_rdata;
        FTensor m_fdata;
        DTensor m_ddata;

        Tensor(RTensor r, FTensor f, DTensor d) : m_rdata{r}, m_fdata{f}, m_ddata{d}
        {}

        friend class Maps;

    public:

        Tensor() : m_rdata{}, m_fdata{}, m_ddata{}
        {}

        explicit Tensor(TKey k) : m_rdata(k), m_fdata(k), m_ddata(k)
        {}

        explicit Tensor(DEG depth) : m_rdata{}, m_fdata{}, m_ddata{}
        {
            using VBP = alg::vectors::VectorBasisProperties<typename RTensor::BASIS>;
            DIMN size = VBP::DegProp::start_of_degree(depth+1);
            m_fdata.maybe_resize(size);
            m_ddata.maybe_resize(size);
        }

        template<typename F>
        void fill_with(F fn, DEG depth)
        {
            using VBP = alg::vectors::VectorBasisProperties<typename RTensor::BASIS>;
            DIMN size = VBP::DegProp::start_of_degree(depth+1);
            std::vector<Sca> tmp;
            tmp.reserve(size);

            typename RTensor::KEY key {};
            auto basis = RTensor::basis;

            for (DIMN i=0; i<size; ++i) {
                Sca val = fn(i);
                m_rdata[key] = val;
                key = basis.nextkey(key);
                tmp.push_back(val);
            }

            REQUIRE CHECK(tmp.size() == size);

            m_fdata.fill_with([&] (DIMN i) { return float(tmp[i]); }, depth);
            m_ddata.fill_with([&] (DIMN i) { return double(tmp[i]); }, depth);
        }

        template <typename F, typename RNG, typename DIST>
        void insert_sparse(F fn, DEG start, DEG end, RNG rng, DIST dist, DIMN max_skip=25)
        {
            using Key = typename RTensor::KEY;
            REQUIRE CHECK(start <= Depth);

            // Construct start key
            Key k;
            for (int i=0; i<start; ++i)
            {
                k.push_back(LET(1));
            }

            auto& basis = RTensor::basis;

            auto advance_key = [&] () {
                DIMN skip = 1 + (rng() % max_skip);
                // skip random number of steps to get the next key
                for (int _i=0; _i<skip; ++_i) { k = basis.nextkey(k); }
            };

            advance_key();
            while (basis.degree(k) <= end) {
                Sca val = dist(rng);

                m_rdata[k] = val;
                m_fdata[k] = float(val);
                m_ddata[k] = double(val);

                advance_key();
            }


        }

//        DECLARE_UMINUS(Tensor)

        DECLARE_BINARY_OPERATOR(Tensor, +, Tensor)
        DECLARE_BINARY_OPERATOR(Tensor, -, Tensor)

        DECLARE_BINARY_OPERATOR(Tensor, *, Sca)
        DECLARE_BINARY_OPERATOR(Tensor, /, Rat)

        DECLARE_BINARY_OPERATOR(Tensor, *, Tensor)

        DECLARE_INPLACE_BINARY_OP(Tensor, +=, Tensor)
        DECLARE_INPLACE_BINARY_OP(Tensor, -=, Tensor)

        DECLARE_INPLACE_BINARY_OP(Tensor, +=, Sca)
        DECLARE_INPLACE_BINARY_OP(Tensor, /=, Rat)

        DECLARE_INPLACE_BINARY_OP(Tensor, *=, Tensor)

        friend Tensor exp(const Tensor& arg)
        {
            return Tensor(exp(arg.m_rdata), exp(arg.m_fdata), exp(arg.m_ddata));
        }

        friend Tensor log(const Tensor& arg)
        {
            return Tensor(log(arg.m_rdata), log(arg.m_fdata), log(arg.m_ddata));
        }

        friend Tensor inv(const Tensor& arg)
        {
            return Tensor(inv(arg.m_rdata), inv(arg.m_fdata), inv(arg.m_ddata));
        }

        friend void check(const Tensor& arg, float float_error=2.0e-4, double double_error=2.0e-12)
        {
            using Key = typename RTensor::KEY;
            DIMN ref_size = arg.m_rdata.size();
            CHECK_EQUAL(arg.m_fdata.size(), ref_size);
            CHECK_EQUAL(arg.m_ddata.size(), ref_size);

            std::vector<std::pair<Key, Sca>> ref_data;
            ref_data.assign(arg.m_rdata.begin(), arg.m_rdata.end());
            std::sort(ref_data.begin(), ref_data.end(),
                      [](const std::pair<Key, Sca> &lhs,
                         const std::pair<Key, Sca> &rhs) -> bool {
                          return lhs.first < rhs.first;
                      }
            );

            auto ref_d = ref_data.begin();


            Key r_key, f_key, d_key;
            Sca r_val;
            float f_rval, f_val;
            double d_rval, d_val;

            for (DIMN i=0; i<ref_size; ++i)
            {
                assert (ref_d != ref_data.end());

                r_key = ref_d->first; r_val = ref_d->second;
                auto f_d = arg.m_fdata.find(r_key);
                auto d_d = arg.m_ddata.find(r_key);

                REQUIRE CHECK(f_d != arg.m_fdata.end());
                REQUIRE CHECK(d_d != arg.m_ddata.end());

                f_key = f_d->first; f_val = f_d->second;
                d_key = d_d->first; d_val = d_d->second;
                f_rval = float(r_val);
                d_rval = double(r_val);

                REQUIRE CHECK_EQUAL(r_key, f_key);
                REQUIRE CHECK_CLOSE(f_rval, f_val, float_error);

                REQUIRE CHECK_EQUAL(r_key, d_key);
                REQUIRE CHECK_CLOSE(d_rval, d_val, double_error);

                ++ref_d;
            }


        }

    };

    class Maps
    {
        typedef typename RationalF::MAPS RMaps;
        typedef typename FloatF::MAPS FMaps;
        typedef typename DoubleF::MAPS DMaps;

        RMaps m_rmaps;
        FMaps m_fmaps;
        DMaps m_dmaps;


    public:

        Tensor l2t(const Lie& arg)
        {
            return Tensor(
                    m_rmaps.l2t(arg.m_rdata),
                    m_fmaps.l2t(arg.m_fdata),
                    m_dmaps.l2t(arg.m_ddata)
            );
        }

    };


#undef DECLARE_BINARY_OPERATOR
#undef DECLARE_INPLACE_BINARY_OP
#undef DECLARE_UMINUS

    AccuracyTest()
    {}

};





#endif //LIBALGEBRAUNITTESTS_ACCURACY_TEST_H
