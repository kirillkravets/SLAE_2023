#ifndef BCG_METHOD
#define BCG_METHOD

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include <fstream>



// BiCG method
template<typename T>
std::vector<T> BCG_method(const CsrMatrix<T> A,  const std::vector<T>& x0, const std::vector<T>& b, T tolerance)  {

    CsrMatrix<T> At = A.transpose();

    std::vector<T> x = x0;

    std::vector<T> r0 = A * x - b;

    T norm = Third_Norm(r0);

    std::vector<T> r = A * x - b;
    std::vector<T> r_wave = r;

    std::vector<T> p = r;
    std::vector<T> p_wave = r;
    
    std::vector<T> z;
    std::vector<T> z_wave;

    T rho;
    T rho_previous;

    T q;
    T theta;

    size_t count = 0;

    while (norm > tolerance) {

        
        for (std::size_t i = 0; i < x.size(); ++i) 
        {

            rho = r * r_wave;

            if (rho == 0) 
                break;

            if( i == 0){
                rho_previous = rho;
            }

            else {
                theta = rho / rho_previous;

                p = r + theta * p;
                p_wave = r_wave + theta * p_wave;
            }
  
            z = A * p;
            z_wave = At * p_wave;

            q = rho / (r_wave * z);

            x = x - q * p;
            r = r - q * z;

            r_wave = r_wave - q * z_wave;

            rho_previous = rho;
            
            norm = Third_Norm(r);
            count++;
        }
    }

    return x;
}

#endif