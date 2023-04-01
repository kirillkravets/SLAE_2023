#ifndef GAUSS_ZEIDEL
#define GAUSS_ZEIDEL

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"

#include<fstream>


template <typename T>
std::vector<T> GaussSeidelMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r){
   
    std::ofstream fout;
    std::ofstream fout1;

    fout.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_gauss_zeidel_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_gauss_zeidel_lnr.txt");
    
    std::vector<T> r_i;
    std::vector<T> x = x0;

    std::vector<T> A_elems = A.Get_Elements();
    std::vector<std::size_t> A_col_ind = A.Get_Col_Ind();
    std::vector<std::size_t> A_amount_of_elems = A.Get_Amount_Of_Elems(); 

    T norm = r + 1;

    std::size_t it = 1;

    while (norm > r)
    {
        fout << it << '\n';

        r_i = b - A * x;

        for(std::size_t i = 0; i < x.size(); i++){
            
            x[i] = b[i];
            
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

        norm = Third_Norm(r_i);

        fout1 << log(norm) << '\n';

        it++;  
    }

    fout.close();
    fout1.close();

    return x;
}

#endif