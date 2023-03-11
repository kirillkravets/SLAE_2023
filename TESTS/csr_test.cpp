#include <gtest/gtest.h>
#include "../source/Matrixes/csr_matrix/csr_matrix.hpp"

#include <iostream>

using namespace struct_DOC;


TEST(test_1, subtest_3x3){
    
 
    vector<DOC<double>> vec_of_matrix({{0, 0, 6},{1,2, 4},{5, 3, 8}});

    CsrMatrix<double> matrix = CsrMatrix(vec_of_matrix);

    for(size_t i = 0; i < matrix.SizeCol(); i++ )
    {
        for(size_t j = 0; j < matrix.SizeStr(); j++)
            std::cout << matrix(i, j) << ' ';
        std::cout << '\n';
    }


    vector<double> v({1, 2, 3, 4});
    
    vector<double> vec = matrix * v;
    for(size_t i = 0; i < matrix.SizeCol(); i++){
        std::cout << vec[i] << ' ';
    }

    
}



int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
