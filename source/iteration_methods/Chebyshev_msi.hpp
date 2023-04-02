#ifndef CHEBYSHEV
#define CHEBYSHEV

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"




template<typename T>
std::vector<T> ChebyshevMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r, T lambda_min, T lambda_max, size_t degree){
    
    std::size_t n = x0.size();
    std::vector<T> z(degree);

    z[0]        = std::cos(M_PI / (2 * n));
    T sin_beta  = std::sin(M_PI / (2 * n));
    T cos_alpha = std::cos(M_PI / n);
    T sin_alpha = std::sin(M_PI / n);
    

    for(std::size_t i = 1; i < degree; i++){
        z[i] = z[i-1] * cos_alpha - sin_beta * sin_alpha;
        sin_beta = z[i] * sin_alpha + cos_alpha * sin_beta;
    }

    for(std::size_t i = 1; i < degree; i++)
    {
        z[i] = 1 / ((lambda_max - lambda_min) / 2 * z[i] + (lambda_max + lambda_min) / 2);
    }

    std::size_t n0 = degree;
    
    std::vector<std::size_t> indexes_z(degree);

    for(std::size_t i = 1; i < degree; i*=2){
        for(std::size_t j = 0; j  < degree ; j += n0)

            indexes_z[j + n0/2] = 2 * i - indexes_z[j] - 1;

        n0 /= 2;    
    }

    std::vector<T> x = x0;
    std::vector<T> r_i(x.size());
    
    T norm = r + 1;
    
    r_i = b - A * x;
    while(norm > r)
    { 

        for(std::size_t i = 0; i < degree; i++){
            
            x = x + z[indexes_z[i]] * r_i;

            r_i = b - A * x;

        }

        norm = Third_Norm(r_i);
    }


    return x;
}


#endif