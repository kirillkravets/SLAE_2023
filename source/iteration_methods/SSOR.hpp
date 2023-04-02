#ifndef SSOR
#define SSOR

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"


template<typename T>
std::vector<T> SSOR_method_fragment(const std::vector<T>& A_elems, const std::vector<std::size_t>& A_col_ind, const std::vector<std::size_t>& A_amount_of_elems, const std::vector<T> x0, const std::vector<T>& b, T omega)
{
    std::vector<T> x = x0;
    double x_i;
        
    for(std::size_t i = 0; i < x.size(); i++)
    {
        x_i = x[i];
        x[i] = omega * b[i];
        
        T diag;

        for(std::size_t j = A_amount_of_elems[i]; j < A_amount_of_elems[i+1]; j++)
        {   
            if(A_col_ind[j] == i)
            {
                diag = A_elems[j]; 
                continue;
            } 

            x[i] -= omega * A_elems[j] * x[A_col_ind[j]];
        }

        x[i] /= diag;
        x[i] += (1 - omega) * x[i];
    }

    for(std::size_t i = x.size() - 2; i >= 0; i--)
    {

        x_i = x[i];

        x[i] = omega * b[i];
        
        T diag;

        for(std::size_t j = A_amount_of_elems[i]; j < A_amount_of_elems[i+1]; j++)
        {   
            if(A_col_ind[j] == i)
            {
                diag = A_elems[j]; 
                continue;
            } 

            x[i] -= omega * A_elems[j] * x[A_col_ind[j]];
        }

        x[i] /= diag;
        x[i] += (1 - omega) * x[i];
    }
    
    return x;
}



#endif