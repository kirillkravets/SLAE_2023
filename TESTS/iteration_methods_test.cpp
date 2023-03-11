#include <gtest/gtest.h>
#include "../source/Matrixes/csr_matrix/csr_matrix.hpp"
#include "../source/iteration_methods/Gauss_Zeidel.hpp"
#include "../source/iteration_methods/Method_Jacobi.hpp"
#include "../source/iteration_methods/Method_simple_iterations.hpp"

#include <iostream>

using namespace struct_DOC;


TEST(test_1, task_3){
    
 
    std::vector<DOC<double>> vec_of_matrix({
        {0, 0, 10 },\
        {0, 1, 1  },\
        {1, 0, 1  },\
        {1, 1, 7 },\
        {2, 1, 0.1},\
        {2, 2, 1  }
    });

    CsrMatrix<double> A = CsrMatrix(vec_of_matrix);

    std::vector<double> x0({0, 0, 0});
    std::vector<double> b({20, 30, 1});

    double r = 1e-12;
    // double t = 1e-4;
    
    std::vector<double> x1 = Method_Yacobi(A, b, x0, r);

    for(size_t i = 0; i < x1.size(); i++){
        std::cout << x1[i] << ' ';
    }

    std::cout << "\n\n";

    std::vector<double> x2 = GaussSeidelMethod(A, b, x0, r);

    for(size_t i = 0; i < x2.size(); i++){
        std::cout << x2[i] << ' ';
    }

    std::cout << "\n\n";


    // std::vector<double> x3 = MethodSimpleIterations(A, b, x0, t, r);
    // for(size_t i = 0; i < x3.size(); i++){
    //     std::cout << x3[i] << ' ';
    // }

    // std::cout << "\n\n";

}



TEST(test_2, task_4){
    
 
    std::vector<DOC<double>> vec_of_matrix({
        {0, 0, 12 },\
        {0, 1, 17 },\
        {0, 2, 3  },\
        {1, 0, 17 },\
        {1, 1, 15825},\
        {1, 2, 28 },\
        {2, 0, 3  },\
        {2, 1, 28 },\
        {2, 2, 238},\
    });

    CsrMatrix<double> A = CsrMatrix(vec_of_matrix);

    std::vector<double> x0({1, 1, 1});
    std::vector<double> b({1, 2, 3});

    std::vector<double> c({1e10, 1e-12, 1e5});
    

    double r = 1e-12;
    double t = 1e-4;
    
    std::vector<double> x1 = Method_Yacobi(A, b, x0, r);

    // for(size_t i = 0; i < x1.size(); i++){
    //     std::cout << x1[i] << ' ';
    // }

    // double min1 = func(A, x1, b);
    // std::cout << "\nf(x_min): "<<  min1 << "\n\n";
    
    std::vector<double> x2 = GaussSeidelMethod(A, b, x0, r);

    // for(size_t i = 0; i < x2.size(); i++){
    //     std::cout << x2[i] << ' ';
    // }

    // double min2 = func(A, x2, b);
    // std::cout << "\nf(x_min): "<<  min2 << "\n\n";


    std::vector<double> x3 = MethodSimpleIterations(A, b, x0, t, r);
    // for(size_t i = 0; i < x3.size(); i++){
    //     std::cout << x3[i] << ' ';
    // }

    // double min3 = func(A, x3, b);
    // std::cout << "\nf(x_min): "<<  min3 << "\n\n";



}



int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
