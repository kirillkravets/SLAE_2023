#ifndef TRIDIAGONAL_MATRIX_HPP
#define TRIDIAGONAL_MATRIX_HPP

#include <vector>
using std::vector;

namespace triple{
    template <typename T>
    struct Triple{
        
        T a;
        T b;
        T c;
    };
};

template <typename T>
class TridiagonalMatrix{
private:

    vector<triple::Triple<T>> triples;

    std::size_t matrix_order;

public:

    TridiagonalMatrix(const vector<triple::Triple<T>>& _triples)
    {
        triples = _triples;
        matrix_order = triples.size();
    }

    std::size_t Get_Order() const
    {
        return matrix_order;
    }

    triple::Triple<T> Get_Triple(std::size_t i) const
    {
        return triples[i];
    }

    triple::Triple<T> operator [] (int i) const{
        return triples[i];
    }
};

#endif /*TRIDIAGONAL_MATRIX_HPP*/