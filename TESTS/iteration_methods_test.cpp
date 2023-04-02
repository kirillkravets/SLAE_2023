#include <gtest/gtest.h>
#include "../source/Matrixes/csr_matrix/csr_matrix.hpp"
#include "../source/iteration_methods/Gauss_Zeidel.hpp"
#include "../source/iteration_methods/Method_Jacobi.hpp"
#include "../source/iteration_methods/Method_simple_iterations.hpp"
#include "../source/iteration_methods/Chebyshev_msi.hpp"
#include "../source/iteration_methods/Gauss_Zeidel_Sym.hpp"

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
    std::vector<double> solution({1.594, 4.058, 0.594});
    
    std::vector<double> x1 = Method_Yacobi(A, b, x0, r);

    for(std::size_t i = 0; i < 3; i++){
        ASSERT_NEAR(x1[i], solution[i], 0.001)
        << ">> FAILED YACOBI TEST COORD: " << i << '!' << std::endl;  
    }


    std::vector<double> x2 = GaussSeidelMethod(A, b, x0, r);

    for(std::size_t i = 0; i < 3; i++){
        ASSERT_NEAR(x2[i], solution[i], 0.001)
        << ">> FAILED GAUSS-SEIDEL TEST COORD: " << i << '!' << std::endl;  
    }
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
    
    std::vector<double> solution({0.080408, 0.000019, 0.011589});

    std::vector<double> x1 = Method_Yacobi(A, b, x0, r);


    for(std::size_t i = 0; i < 3; i++){
        ASSERT_NEAR(x1[i], solution[i], 1e-6)
        << ">> FAILED YACOBI TEST COORD: " << i << '!' << std::endl;  
    }

    std::vector<double> x2 = GaussSeidelMethod(A, b, x0, r);

    for(std::size_t i = 0; i < 3; i++){
        ASSERT_NEAR(x2[i], solution[i], 1e-6)
        << ">> FAILED GAUSS-SEIDEL TEST COORD: " << i << '!' << std::endl;  
    }

    std::vector<double> x3 = MethodSimpleIterations(A, b, x0, t, r);
    

    for(std::size_t i = 0; i < 3; i++){
        ASSERT_NEAR(x3[i], solution[i], 1e-6)
        << ">> FAILED SIMPLE ITERATIONS TEST COORD: " << i << '!' << std::endl;  
    }

}

TEST(test_3, task_0){
    std::vector<DOC<double>> vec_of_matrix({
        {0, 0, 10 },\
        {0, 1, -0.5 },\
        {1, 0, -0.5 },\
        {1, 1, 10}
    });

    CsrMatrix<double> A = CsrMatrix(vec_of_matrix);

    std::vector<double> x0({0, 0});
    std::vector<double> b({1, 3});

    std::vector<double> c({1e10, 1e-12, 1e5});

    std::vector<double> solution({0.115288, 0.305764});
    double r = 1e-12;

    std::vector<double> x3 = ChebyshevMethod(A, b, x0, r, 9.5, 10.5, 32);

    std::cout << std::endl;


    for(std::size_t i = 0; i < 2; i++){
        ASSERT_NEAR(x3[i], solution[i], 1e-3)
        << ">> FAILED SIMPLE ITERATIONS TEST COORD: " << i << '!' << std::endl;  
    }

    std::vector<double> x4 = GaussSeidelSymMethod(A, b, x0, r, 10.5 - 9.5);
    for(std::size_t i = 0; i < x4.size(); i++){
        ASSERT_NEAR(x4[i], solution[i], 1e-3)
        << ">> FAILED SIMPLE ITERATIONS TEST COORD: " << i << '!' << std::endl;  
    }
}



int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
