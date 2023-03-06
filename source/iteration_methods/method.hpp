#include "../csr_matrix/csr_matrix.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include <cmath>
using std::size_t;

vector<double> operator+(const vector<double>& a, const vector<double>& b)  {
    
    vector<double> result;
    result.reserve(a.size());

    for (size_t i = 0; i < a.size(); i++){
        result.push_back(a[i] + b[i]);
    }

    return result;
}

vector<double> operator-(const vector<double>& a, const vector<double>& b)  {
    
    vector<double> result;
    result.reserve(a.size());

    for (size_t i = 0; i < a.size(); i++){
        result.push_back(a[i] - b[i]);

    }

    return result;
}

vector<double> operator*(const vector<double>& a, const double digit){
    
    vector<double> result;
    result.reserve(a.size());

    for (size_t i = 0; i < a.size(); i++){
        result.push_back(a[i] * digit); 
    }

    return result;
}

vector<double> operator*(const double digit, const vector<double>& a){
    
    return a * digit;
}

double Norm_Of_Vec(const vector<double>& vec){

    double result = 0;


    for(size_t i = 0; i < vec.size(); i++){
        
        result += vec[i] * vec[i];
    }

    return result;
}



vector<double> MethodSimpleIterations(const CsrMatrix<double>& A, const vector<double>& b, const vector<double>& x0, double t, double r) {
    std::ofstream fout;
    std::ofstream fout1;
    fout.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_msi_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_msi_lnr.txt");

    vector<double> x = x0;
    vector<double> r_i;
    
    double norm = r + 1;
    
    size_t it = 1;

    while(norm > r)
    { 
        fout << it << '\n';
        r_i = b - A * x;
        
        x = x + t * r_i;

        norm = Norm_Of_Vec(r_i);
        fout1 << 0.5 * log(norm) << '\n';
        //std::cout << norm << '\n';
        it++;
    }
    fout.close();
    fout1.close();
    std::cout << "SimpleIterations complete!\n";


    return x;
}


vector<double> YacobiMethod(const CsrMatrix<double>& A, const vector<double>& b, const vector<double>& x0, double r){
    std::ofstream fout;
    std::ofstream fout1;
    fout.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_jacobi_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_jacobi_lnr.txt");

    std::size_t size_vectors = x0.size();
    
    vector<double> r_i;
    vector<double> x = x0;
    
    double norm = r + 1;
    size_t it = 1;
    while (norm > r)
    {
        fout << it << '\n';
        r_i = b - A * x;

        for(size_t i = 0; i < size_vectors; i++){
            
            double x_i = b[i];
            
            for(size_t j = 0; j < size_vectors; j++)
            {
                if(i != j){
                    x_i -= A(i, j) * x[j];
                }
            }

            x[i] = x_i / A(i, i);
        }
        
        norm = Norm_Of_Vec(r_i);
        fout1 << 0.5 * log(norm) << '\n';

        it++;
    }

    fout.close();
    fout1.close();

    std::cout << "Jacobi complete!\n";
    return x;
}


vector<double> GaussSeidelMethod(const CsrMatrix<double>& A, const vector<double>& b, const vector<double>& x0, double r){
   
    std::ofstream fout;
    std::ofstream fout1;
    fout.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_gauss_zeidel_it.txt");
    fout1.open("/home/kirill/vs codes c++/SLAE_projects/files_txt/tsk4_gauss_zeidel_lnr.txt");


    std::size_t size_vectors = x0.size();
    
    vector<double> r_i;
    vector<double> x = x0;
    vector<double> x1;

    double norm = r + 1;
    size_t it = 1;

    while (norm > r)
    {
        fout << it << '\n';

        x1 = x;
        r_i = b - A * x1;

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
            x[i] = x1[i];

        }

        norm = Norm_Of_Vec(r_i);

        fout1 << 0.5 * log(norm) << '\n';

        it++;
        
        
        
    }

    fout.close();
    fout1.close();

    std::cout << "Gauss-Seidel complete!\n";

    return x1;
}
 
double func(const CsrMatrix<double>& A, const vector<double>& x, const vector<double>& b)
{
    double result = 0;

    for(size_t i = 0; i < b.size(); i++){
        for(size_t j = 0; j < b.size(); j++)
        {
            result += x[i] * A(i, j) * x[j];
        }

        result += x[i] * b[i];
    }

    return result;
}