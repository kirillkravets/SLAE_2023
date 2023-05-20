#ifndef BCG_METHOD
#define BCG_METHOD

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include <fstream>




template<typename T>
std::vector<T> BCG_method(const CsrMatrix<T>& A, const std::vector<T>& x_0, const std::vector<T>& b, T r0)
{
    
    std::vector<T> x = x_0;
   
    std::vector<T> r = A * x - b;
    
    T norm = Third_Norm(r);

    CsrMatrix<T> A_transp = A.transpose();

    std::size_t count = 0;
   
    while (norm > r0) {
        
        std::vector<T> r1 = A * x - b;
        std::vector<T> alpha_v = r1;

        std::vector<T> rr = r1;
        std::vector<T> RR;

        std::vector<T> beta_v = r1;

        T q, t, next_r0, previous_r0;
        
        
        for (int i = 0; i < A.SizeStr(); ++i) {
            
            next_r0 = r1 * alpha_v;
          
            if (next_r0 == 0){ 
                break;
            }
            
            if (i != 0) {
                t = next_r0 / previous_r0;
                rr = r1 + t * rr;
                beta_v = alpha_v + t * beta_v;
            }
           
            else{
                previous_r0 = next_r0;
            }


            RR = A * rr;

            q = next_r0 / (alpha_v * RR);
            previous_r0 = next_r0;
          
            alpha_v = alpha_v - q * (A_transp * beta_v);
       
            x = x - q * rr;
            
            r1 = r1 - q * RR;            

            norm = Third_Norm(r);
            count++;
        }
    }
    return x;
}


#endif