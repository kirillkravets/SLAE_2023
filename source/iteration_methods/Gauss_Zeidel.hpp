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


    std::size_t size_vectors = x0.size();
    
    std::vector<T> r_i;
    std::vector<T> x = x0;
    std::vector<T> x1;
    x1.resize(x.size());

    T norm = r + 1;

    std::size_t it = 1;

    while (norm > r)
    {
        fout << it << '\n';

        r_i = b - A * x;

        for(std::size_t i = 0; i < size_vectors; i++){
            
            x1[i] = b[i];
            
            for(std::size_t j = i + 1; j < size_vectors; j++)
            {               
                x1[i] -= A(i, j) * x[j];
            }

            for(std::size_t j = 0; j < i; j ++)
            {
                x1[i] -= A(i, j) * x1[j];
            }

            x1[i] /= A(i, i);

        }

        x = x1;

        norm = Third_Norm(r_i);

        fout1 << log(norm) << '\n';

        it++;  
    }

    fout.close();
    fout1.close();

    return x;
}

#endif