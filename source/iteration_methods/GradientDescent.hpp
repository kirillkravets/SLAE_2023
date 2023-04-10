#ifndef GRADIENT_DESCENT_METHOD
#define GRADIENT_DESCENT_METHOD

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include <fstream>

template <typename T>
std::vector<T> GradientDescentMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T  r, T tau, std::string txt1, std::string txt2) {
    
    std::ofstream fout;
    std::ofstream fout1;

    std::string dir_it = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/" + txt1 + ".txt";
    std::string dit_lnr = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/" + txt2 + ".txt";
    
    std::string dir_x1 = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/proection_files/" + txt1 + "_x1.txt";
    std::string dir_x4 = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/proection_files/" + txt2 + "_x4.txt";
    

    fout.open(dir_it);
    fout1.open(dit_lnr);

    std::ofstream foutx1;
    std::ofstream foutx4;

    foutx1.open(dir_x1);
    foutx4.open(dir_x4);

    std::size_t it = 1;


    std::vector<T> x = x0;
    std::vector<T> r_i(b.size());

    r_i = b - A * x;

    T norm = Third_Norm(r_i);
    
    while(norm > r)
    {   
        x = x + r_i * tau;

        r_i = b - A * x;

        norm = Third_Norm(r_i);

        fout << it << '\n';
        fout1 << log(norm) << '\n';
        it++;

        foutx1 << x[0] << '\n';
        foutx4 << x[x.size() - 1] << '\n';
    }

    foutx1 << x[0] << '\n';
    foutx4 << x[x.size() - 1] << '\n';

    fout.close();
    fout1.close();

    foutx1.close();
    foutx4.close();

    return x;
}


#endif