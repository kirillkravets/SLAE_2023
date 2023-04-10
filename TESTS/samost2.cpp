#include <gtest/gtest.h>
#include "../source/Matrixes/csr_matrix/csr_matrix.hpp"
#include "../source/iteration_methods/Gauss_Zeidel.hpp"
#include "../source/iteration_methods/Method_Jacobi.hpp"
#include "../source/iteration_methods/Method_simple_iterations.hpp"
#include "../source/iteration_methods/Chebyshev_msi.hpp"
#include "../source/iteration_methods/Gauss_Zeidel_Sym.hpp"
#include "../source/iteration_methods/SSOR.hpp"
#include "../source/iteration_methods/SteepestDescentMethod.hpp"
#include "../source/iteration_methods/GradientDescent.hpp"

#include <iostream>

using namespace struct_DOC;


TEST(test_1, TASK_1){
    
    //double a = 6, b = 17, c = 1;

    std::vector<DOC<double>> vec_of_matrix;
    vec_of_matrix.reserve(289 + (288  + (289 - 17)) * 2);

    for(std::size_t i = 0; i < 289; i++){
        
        if( i >= 17)
            vec_of_matrix.push_back({i, i - 17, 6});

        if( i  >=  1)
            vec_of_matrix.push_back({i, i - 1, 6});
        
        vec_of_matrix.push_back({i, i, 17 * 2});
        
        if(i + 1 <= 288)
            vec_of_matrix.push_back({i, i + 1, 6});

        if(i + 17 <= 288)
        vec_of_matrix.push_back({i, i + 17, 6});
    }

    CsrMatrix<double> A = CsrMatrix(vec_of_matrix);

    std::vector<double> x0(289, 0);

    std::vector<double>  b(289, 1);


    double r = 1e-13;

    double lambda_max = 2 * (17 + 2 * 6 * cos(M_PI / 18));
    double lambda_min = 2 * (17 - 2 * 6 * cos(M_PI / 18));

    std::vector<double> x = SteepestDescentMethod(A, x0, b, r);

    double tau1 = 1/ lambda_max;
    double tau2 = 2/ (lambda_max + lambda_min);

    
    std::vector<double> x1 = GradientDescentMethod(A, b, x0, r, tau1, "1_it", "1_ln");

    std::vector<double> x2 = GradientDescentMethod(A, b, x0, r, tau2, "2_it", "2_ln");

    std::vector<double> x3 = ChebyshevMsiMethod(A, b, x0, r, lambda_min, lambda_max, pow(2, 8), "3_it", "3_ln");

    std::vector<double> x4 = SSOR_method(A, b, x0, r, tau2, 1.2, "4_it", "4_ln");

    // for(std::size_t i = 0; i < x1.size(); i++)
    // {
    //     std::cout << x1[i] << ' ' << x2[i] << ' ' << x3[i] << ' ' << x4[i] << '\n';
    // } 

    std::ofstream fout;
    std::ofstream fout1;

    std::string dir_lambda = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/lambda.txt";
    std::string dir_time = "/home/kirill/vs codes c++/tyyh5etyh5trh/SLAE_2023/TESTS/test_08_04/time.txt";
    
    fout.open(dir_lambda);
    fout1.open(dir_time);

    lambda_min += lambda_min * static_cast<double>(50) / 1000;
    lambda_max -= lambda_max * static_cast<double>(50) / 1000;

    for(int i = 0; i < 94; i++){
        std::string lambda_name = "cheb" + std::to_string(i) + "_lambda";
        std::string time_name   = "cheb" + std::to_string(i) + "_time";

        lambda_min -= lambda_min * static_cast<double>(i) / 1000;
        lambda_max += lambda_max * static_cast<double>(i) / 1000;

        double time = ChebyshevMsiTimeMethod(A, b, x0, r, lambda_min, lambda_max, pow(2, 8)); 

        fout << (lambda_max - lambda_min) << '\n';
        fout1 << time << '\n';

        std::cout << i << ":    " << time <<  "     " <<lambda_max << "  " << lambda_min << '\n';
    }


    fout.close();
    fout1.close();
}




TEST(test_2, TASK_2){
    
    //double a = 6, b = 17, c = 1;
 
    std::vector<DOC<double>> vec_of_matrix({
        {0, 0, 6 },\
        {1, 1, 7  },\
        {2, 2, 8  },\
        {3, 3, 9 },\
    });

    CsrMatrix<double> A = CsrMatrix(vec_of_matrix);

    std::vector<double> x0({0, 0, 0, 0});
    std::vector<double> b({5, 5, 5, 5});
    std::vector<double> c({6, 6, 6, 6});
    
    double r = 1e-13;

    std::vector<double> x = SteepestDescentMethod(A, x0, b, r);
    
    for(std::size_t i = 0; i < x0.size(); i++){
        std::cout << x[i] << ' ';
    }

    std::cout << '\n';

    double f_min = func(A, x, b) + 6;

    std::cout << f_min;

    std::cout << '\n';
    
}




int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
