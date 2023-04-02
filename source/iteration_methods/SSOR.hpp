#ifndef SSOR
#define SSOR

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"


template<typename T>
std::vector<T> SSOR_method_fragment(const std::vector<T>& A_elems, const std::vector<std::size_t>& A_col_ind, const std::vector<std::size_t>& A_amount_of_elems, const std::vector<T>& x0, const std::vector<T>& b, T omega)
{
    std::vector<T> x = x0;
    
    T x_i;
        
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
        x[i] += (1 - omega) * x_i;
    }

    for(std::size_t i = x.size() - 2; i > 0; i--)
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

            x[i] -= omega * (A_elems[j] * x[A_col_ind[j]]);
        }

        x[i] /= diag;
        x[i] += (1 - omega) * x_i;

    }
    
    return x;
}



template <typename T>
std::vector<T> SSOR_method(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r, T rho, T omega){
   
    std::vector<T> r_i;
    
    std::vector<T>           A_elems           = A.Get_Elements();
    std::vector<std::size_t> A_col_ind         = A.Get_Col_Ind();
    std::vector<std::size_t> A_amount_of_elems = A.Get_Amount_Of_Elems(); 

    std::vector<T> y0 = x0;
    std::vector<T> y = SSOR_method_fragment(A_elems, A_col_ind, A_amount_of_elems, y0, b, omega);

    std::vector<T> gzsf = SSOR_method_fragment(A_elems, A_col_ind, A_amount_of_elems, y, b, omega);

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

        // for(std::size_t i = 0; i < y0.size(); i++){
        //     std::cout << y0[i] << ' ';
        // }
        // std::cout << '\n';

        r_i = b - A * y;

        norm = Third_Norm(r_i);

        gzsf = SSOR_method_fragment(A_elems, A_col_ind, A_amount_of_elems, y0, b, omega);
    }

    return y;
}


#endif