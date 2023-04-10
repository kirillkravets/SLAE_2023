#ifndef STEEPESTDESCENTMETHOD
#define STEEPESTDESCENTMETHOD

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include <fstream>
#include <cmath>

template<typename T>
std::vector<T> SteepestDescentMethod(const CsrMatrix<T>& A, const std::vector<T>& x_0, const std::vector<T>& b, T r, std::string txt1, std::string txt2)
{
    std::string dir_x1 = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/proection_files/" + txt1 + "_x1.txt";
    std::string dir_x4 = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/proection_files/" + txt2 + "_x4.txt";

    std::ofstream foutx1;
    std::ofstream foutx4;

    foutx1.open(dir_x1);
    foutx4.open(dir_x4);

    T alpha;

    std::vector<T> x = x_0;

    std::vector<T> r_i(b.size());
    
    r_i = b - A * x;

    std::vector<T> A_r_i(b.size());

    T norm = Third_Norm(r_i);

    while (norm > r) {

        A_r_i = A * r_i;

        alpha = (r_i * r_i) / (r_i * A_r_i);
        
        x = x + alpha * r_i;

        r_i = r_i - alpha * A_r_i;
        
        norm = Third_Norm(r_i);

        foutx1 << x[0] << '\n';
        foutx4 << x[x.size() - 1] << '\n';
        
    }

    foutx1 << x[0] << '\n';
    foutx4 << x[x.size() - 1] << '\n';

    foutx1.close();
    foutx4.close();

    return x;
}


#endif