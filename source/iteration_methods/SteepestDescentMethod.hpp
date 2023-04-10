#ifndef STEEPESTDESCENTMETHOD
#define STEEPESTDESCENTMETHOD

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"

#include <cmath>

template<typename T>
std::vector<T> SteepestDescentMethod(const CsrMatrix<T>& A, const std::vector<T>& x_0, const std::vector<T>& b, T r)
{
    T alpha;

    std::vector<T> x = x_0;

    std::vector<T> r_i(b.size());
    
    r_i = b - A * x;

    std::vector<T> A_r_i(b.size());

    T norm = Third_Norm(r_i);

    while (norm > r) {

        A_r_i = A * r_i;

        alpha = (r_i * r_i) / (r_i * A_r_i);
        
        x = x + alpha * r_i;

        r_i = r_i - alpha * A_r_i;
        
        norm = Third_Norm(r_i);
    }
    return x;
}



// template <typename T>
// T valuefunc(const CsrMatrix<T>& A, const std::vector<T>& x, const std::vector<T>& b)
// {
//     T result = 0;

//     for(std::size_t i = 0; i < b.size(); i++){
//         for(std::size_t j = 0; j < b.size(); j++)
//         {
//             result += x[i] * A(i, j) * x[j];
//         }

//         result += x[i] * b[i];
//     }

//     return result;
// }


#endif