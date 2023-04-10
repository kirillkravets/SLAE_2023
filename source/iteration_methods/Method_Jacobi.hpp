#ifndef METHOD_JACOBI
#define METHOD_JACOBI

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"

#include<fstream>


template <typename T>
std::vector<T> Method_Yacobi(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, const T r){
    
    std::ofstream fout2;
    std::ofstream fout3;
    
    fout2.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_jacobi_it.txt");
    fout3.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_jacobi_lnr.txt.txt");

    std::size_t size_vectors = x0.size();
    
    std::vector<T> r_i;
    std::vector<T> x = x0;
    
    T norm = r + 1;
    std::size_t it = 1;
    
    while (norm > r)
    {
        fout2 << it << '\n';
        r_i = b - A * x;

        for(std::size_t i = 0; i < size_vectors; i++){
            
            T x_i = b[i];
            
            for(std::size_t j = 0; j < size_vectors; j++)
            {
                if(i != j){
                    x_i -= A(i, j) * x[j];
                }
            }

            x[i] = x_i / A(i, i);
        }
        
        norm = Third_Norm(r_i);
        fout3 << log(norm) << '\n';

        it++;
    }

    fout2.close();
    fout3.close();

    return x;
}



#endif