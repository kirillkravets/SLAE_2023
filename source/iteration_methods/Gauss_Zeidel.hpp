#include "../csr_matrix/csr_matrix.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include <cmath>
using std::size_t;


vector<double> GaussSeidelMethod(const CsrMatrix<double>& A, const vector<double>& b, const vector<double>& x0, double r){
   
    std::ofstream fout;
    std::ofstream fout1;
    fout.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_gauss_zeidel_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/SLAE_2023/SLAE_2023/files_txt/tsk4_gauss_zeidel_lnr.txt");


    std::size_t size_vectors = x0.size();
    
    vector<double> r_i;
    vector<double> x = x0;
    vector<double> x1;
    x1.resize(x.size());

    double norm = r + 1;
    size_t it = 1;

    while (norm > r)
    {
        fout << it << '\n';

        r_i = b - A * x;

        for(size_t i = 0; i < size_vectors; i++){
            
            x1[i] = b[i];
            
            for(size_t j = i + 1; j < size_vectors; j++)
            {               
                x1[i] -= A(i, j) * x[j];

            }

            for(size_t j = 0; j < i; j ++){
                x1[i] -= A(i, j) * x1[j];
            }

            x1[i] /= A(i, i);

        }

        x = x1;

        norm = Norm_Of_Vec(r_i);

        fout1 << 0.5 * log(norm) << '\n';

        it++;  
    }

    fout.close();
    fout1.close();

    std::cout << "Gauss-Seidel complete!\n";

    return x1;
}
