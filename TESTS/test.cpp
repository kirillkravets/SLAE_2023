#include <gtest/gtest.h>
#include "../source/tridiagonal_matrix/run_through_method.hpp"


TEST(test_1, subtest_3x3){
    
 
    vector<triple::Triple<double>> vec{{0., 22., 7}, {1., -15., 9.}, {122., -61., 0.}};

    vector<double> right_column{30., 20., 5.};

    TridiagonalMatrix<double> matrix_3x3(vec);

    vector<double> result = Run_Through_method(matrix_3x3, right_column);    

    vector<double> real_result{18.795, 0.647, 1213};

    for(std::size_t i = 0; i < 3; i++){
        std::cout << '\n' << result[i] << '\n';
        //ASSERT_NEAR(result[i], real_result[i], 1.0)
        //<< "!!! TEST FAILED ON COORDINATE NUMBER " << i << " !!!" << std::endl;  
    }
    
}



// TEST(test_2, subtest_3x3){
    
//     vector<triple::Triple<float>> vec{{0,3,10}, {-1, 3, -1}, {-15, 2, 0}};

//     vector<float> right_column{-5, 2, 9};

//     TridiagonalMatrix<float> matrix_3x3(vec);

//     vector<float> result = Run_Through_method(matrix_3x3, right_column);    


//     vector<float> real_result{-7 , 12, -21};

//     for(std::size_t i = 0; i < 3; i++){
//         ASSERT_NEAR(result[i], real_result[i], 1e-5)
//         << ">> FAILED TEST COORD: " << i << '!' << std::endl;  
//     }
// }

int main(int argc, char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
