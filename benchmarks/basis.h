

#include <libalgebra/libalgebra.h>

#ifndef USE_FLATMAP
#define MAP_T std::map<KEY, SCA>
#else
#include <boost/container/flat_map.hpp>
#define MAP_T boost::container::flat_map<KEY, SCA>
#endif
typedef alg::DEG DEG;
typedef alg::DIMN DIMN;

template<DIMN D>
class Basis {
public:
    typedef double RATIONAL;
    typedef double SCALAR;
    typedef size_t KEY;
    typedef std::map<size_t, double> MAP;

    static const DEG MAX_DEGREE = 1 + D / 1000;

    // Default constructor
    Basis() {}

    KEY begin() const 
    {
        return 0;
    }

    KEY nextkey(const KEY& k) const
    {
        return k + 1;
    }

    unsigned degree(const KEY &k) {
        if (k == 0) return 0;
        return 1 + k / 1000;
    }

    KEY end() const
    {
        return std::numeric_limits<KEY>::max();
    }

    static constexpr DEG max_dimension()
    {
        return D;
    }

    static DIMN start_of_degree(const DEG& d)
    {
        if (d == 0) return 0;
        return (d - 1) * 1000;
    }

    inline static KEY key_of_index(const DIMN& i)
    {
        return KEY{i};
    }

    inline static bool comp(const KEY& k1, const KEY& k2)
    {
        return k1 <= k2;
    }

    inline static DIMN index_of_key(const KEY& k)
    {
        return k;
    }

    friend std::ostream& operator<<(
        std::ostream &os,
        const std::pair<Basis*, KEY> &t
    ) {
        return os << t.second;
    }

};
