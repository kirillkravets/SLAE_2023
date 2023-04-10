#ifndef METHOD_SIMPLE_ITERATIONS
#define METHOD_SIMPLE_ITERATIONS

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"

#include<fstream>


template <typename T>
std::vector<T> MethodSimpleIterations(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T t, T  r) {
    
    std::ofstream fout;
    std::ofstream fout1;
    
    fout.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_msi_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/source/iteration_methods/iteration_txt_files/tsk4_msi_lnr.txt");

    std::vector<T> x = x0;
    std::vector<T> r_i;
    
    T norm = r + 1;
    
    std::size_t it = 1;

    while(norm > r)
    { 
        fout << it << '\n';
        r_i = b - A * x;
        
        x = x + t * r_i;

        norm = Third_Norm(r_i);
        fout1 << log(norm) << '\n';
        it++;
    }
    
    fout.close();
    fout1.close();

    return x;
}


template <typename T>
T func(const CsrMatrix<T>& A, const std::vector<T>& x, const std::vector<T>& b)
{
    T result = 0;

    for(std::size_t i = 0; i < b.size(); i++){
        for(std::size_t j = 0; j < b.size(); j++)
        {
            result += x[i] * A(i, j) * x[j];
        }

        result += x[i] * b[i];
    }

    return result;
}

#endif