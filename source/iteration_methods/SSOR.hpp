#ifndef SSOR
#define SSOR

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"


template<typename T>
std::vector<T> SSOR_method(const CsrMatrix<T>& A, const std::vector<T> x0, const std::vector<T>& b, T omega)
{
    std::vector<T> x = x0;

    std::vector<T> A_elems = A.Get_Elements();
    std::vector<std::size_t> A_col_ind = A.Get_Col_Ind();
    std::vector<std::size_t> A_amount_of_elems = A.Get_Amount_Of_Elems(); 

        
    for(std::size_t i = 0; i < x.size(); i++)
    {
        
        x[i] = omega * b[i];
        
        T diag;

        for(std::size_t j = A_amount_of_elems[i]; j < A_amount_of_elems[i+1]; j++)
        {   
            if(A_col_ind[j] == i)
            {
                diag = A(i, i); 
                continue;
            } 

            x[i] -= A_elems[j] * x[A_col_ind[j]];
        }

        x[i] /= diag;
    }

    return x;

}

#endif