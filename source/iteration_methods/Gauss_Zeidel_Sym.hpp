#ifndef GAUSS_ZEIDEL_SYM
#define GAUSS_ZEIDEL_SYM

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"

template<typename T>
std::vector<T> GaussZeidelSymFragment(const std::vector<T>& A_elems, const std::vector<std::size_t>& A_col_ind, const std::vector<std::size_t>& A_amount_of_elems, const std::vector<T>& x0, const std::vector<T>& b)
{
    std::vector<T> x = x0;
        
    for(std::size_t i = 0; i < x.size(); i++)
    {
        x[i] = b[i];
        
        T diag;

        for(std::size_t j = A_amount_of_elems[i]; j < A_amount_of_elems[i+1]; j++)
        {   
            if(A_col_ind[j] == i)
            {
                diag = A_elems[j]; 
                continue;
            } 

            x[i] -= A_elems[j] * x[A_col_ind[j]];
        }

        x[i] /= diag;
    }

    for(std::size_t i = x.size() - 2; i > 0; i--)
    {

        x[i] = b[i];
        
        T diag;

        for(std::size_t j = A_amount_of_elems[i]; j < A_amount_of_elems[i + 1]; j++)
        {   
            if(A_col_ind[j] == i)
            {
                diag = A_elems[j]; 
                continue;
            } 

            x[i] -= A_elems[j] * x[A_col_ind[j]];
        }

        x[i] /= diag;
    }
    
    return x;
}


template <typename T>
std::vector<T> GaussSeidelSymMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r, T rho){
   
    std::vector<T> r_i;
    
    std::vector<T>           A_elems           = A.Get_Elements();
    std::vector<std::size_t> A_col_ind         = A.Get_Col_Ind();
    std::vector<std::size_t> A_amount_of_elems = A.Get_Amount_Of_Elems(); 

    std::vector<T> y0 = x0;
    std::vector<T> y = GaussZeidelSymFragment(A_elems, A_col_ind, A_amount_of_elems, y0, b);

    std::vector<T> gzsf = GaussZeidelSymFragment(A_elems, A_col_ind, A_amount_of_elems, y, b);

    T mu0 = 1;
    T mu = 1/rho;

    T norm = r + 1;

    while (norm > r)
    {
        y0 = -mu0 * y + 2 * (mu / rho) * gzsf;

        mu0 = (2 / rho) * mu - mu0; // mu0 = mu_(i-1) -> mu_(i+1)

        y0 = 1/ mu0 * y0;

        y = y0;
        mu = mu0;

        r_i = b - A * y;

        norm = Third_Norm(r_i);

        gzsf = GaussZeidelSymFragment(A_elems, A_col_ind, A_amount_of_elems, y0, b);
    }

    return y;
}



#endif