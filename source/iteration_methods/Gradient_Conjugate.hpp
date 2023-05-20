#ifndef GRADIENT_CONJUGATE
#define GRADIENT_CONJUGATE

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include <fstream>



template<typename T>
std::vector<T> GradientC_method(const CsrMatrix<T>& A, const std::vector<T>& x_0, const std::vector<T>& b, T r0) 
{
    std::vector<T> x = x_0;
    
    std::vector<T> r = A * x - b;
    std::vector<T> p = r;
    
    T alpha, beta;
 
    std::vector<T> P(x_0.size());
 
    std::size_t count = 0;
  
    T norm = Third_Norm(r);

    while (norm > r0) {
        
        for (int i = 0; i < x_0.size(); ++i) 
        {   
            if (!isZero(r)) 
            {    
                P = A * p;
                alpha = r * r / (p * P);
                x = x - alpha * p;
               
                beta = 1 / (r * r);
              
                r = r - alpha * P;
                beta *= (r * r);
              
                p = r + beta * p;
            
                norm = Third_Norm(r);
                count++;
            }

            else
            {
                return x
            }
        }
    }
    return x;
}

#endif
