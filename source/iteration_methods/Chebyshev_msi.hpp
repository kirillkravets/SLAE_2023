#ifndef CHEBYSHEV
#define CHEBYSHEV

#include "../Matrixes/csr_matrix/csr_matrix.hpp"
#include "../tools/vector_overloads.hpp"
#include "../tools/vec_norm.hpp"
#include<fstream>
#include<time.h>


template<typename T>
std::vector<T> ChebyshevMsiMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r, T lambda_min, T lambda_max, std::size_t degree, std::string it_name_txt, std::string lnr_name_txt){
    

    std::ofstream fout;
    std::ofstream fout1;

    std::string dir_it = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/" + it_name_txt + ".txt";
    std::string dit_lnr = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/" + lnr_name_txt + ".txt";
    
    fout.open(dir_it);
    fout1.open(dit_lnr);
    std::size_t it = 1;

    std::size_t n = x0.size();
    std::vector<T> z(degree);

    z[0]        = std::cos(M_PI / (2 * n));
    T sin_beta  = std::sin(M_PI / (2 * n));
    T cos_alpha = std::cos(M_PI / n);
    T sin_alpha = std::sin(M_PI / n);
    

    for(std::size_t i = 1; i < degree; i++){
        z[i] = z[i-1] * cos_alpha - sin_beta * sin_alpha;
        sin_beta = z[i] * sin_alpha + cos_alpha * sin_beta;
    }

    for(std::size_t i = 1; i < degree; i++)
    {
        z[i] = 1 / ((lambda_max - lambda_min) / 2 * z[i] + (lambda_max + lambda_min) / 2);
    }

    std::size_t n0 = degree;
    
    std::vector<std::size_t> indexes_z(degree);

    for(std::size_t i = 1; i < degree; i*=2){
        for(std::size_t j = 0; j  < degree ; j += n0)

            indexes_z[j + n0/2] = 2 * i - indexes_z[j] - 1;

        n0 /= 2;    
    }

    std::vector<T> x = x0;
    std::vector<T> r_i(x.size());
    
    T norm = r + 1;
    
    r_i = b - A * x;

    while(norm > r)
    { 

        for(std::size_t i = 0; i < degree; i++){
            
            x = x + z[indexes_z[i]] * r_i;

            r_i = b - A * x;

            norm = Third_Norm(r_i);

            fout << it << '\n';
            fout1 << log(norm) << '\n';
            it++;
        }


    }

    fout.close();
    fout1.close();

    return x;
}

template<typename T>
double ChebyshevMsiTimeMethod(const CsrMatrix<T>& A, const std::vector<T>& b, const std::vector<T>& x0, T r, T lambda_min, T lambda_max, std::size_t degree){

    clock_t tStart = clock();
    std::size_t n = x0.size();
    std::vector<T> z(degree);

    z[0]        = std::cos(M_PI / (2 * n));
    T sin_beta  = std::sin(M_PI / (2 * n));
    T cos_alpha = std::cos(M_PI / n);
    T sin_alpha = std::sin(M_PI / n);
    

    for(std::size_t i = 1; i < degree; i++){
        z[i] = z[i-1] * cos_alpha - sin_beta * sin_alpha;
        sin_beta = z[i] * sin_alpha + cos_alpha * sin_beta;
    }

    for(std::size_t i = 1; i < degree; i++)
    {
        z[i] = 1 / ((lambda_max - lambda_min) / 2 * z[i] + (lambda_max + lambda_min) / 2);
    }

    std::size_t n0 = degree;
    
    std::vector<std::size_t> indexes_z(degree);

    for(std::size_t i = 1; i < degree; i*=2){
        for(std::size_t j = 0; j  < degree ; j += n0)

            indexes_z[j + n0/2] = 2 * i - indexes_z[j] - 1;

        n0 /= 2;    
    }

    std::vector<T> x = x0;
    std::vector<T> r_i(x.size());
    
    T norm = r + 1;
    
    r_i = b - A * x;

    while(norm > r)
    { 

        for(std::size_t i = 0; i < degree; i++){
            
            x = x + z[indexes_z[i]] * r_i;

            r_i = b - A * x;

            norm = Third_Norm(r_i);
        }


    }

    return (double)(clock() - tStart)/CLOCKS_PER_SEC;
}

#endif