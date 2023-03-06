#ifndef SWEEP_METHOD_HPP
#define SWEEP_METHOD_HPP

#include "tridiagonal_matrix.hpp"
#include <iostream>
#include <cmath>

template <typename T>
vector<T> Run_Through_method(const TridiagonalMatrix<T>& matrix, const vector<T>& right_column)
{
    
    std::size_t order = matrix.Get_Order();

    vector<T> alpha;
    vector<T> beta;

    alpha.resize(order);
    beta.resize(order);
   
    alpha[0] = -matrix[0].c / matrix[0].b;
    beta[0] = right_column[0]/ matrix[0].b;
    

    for(std::size_t i = 1; i < order; i++){
        
        T denominator = matrix[i].a * alpha[i-1] + matrix[i].b;

        alpha[i] = -matrix[i].c / denominator;
        beta[i] = (right_column[i] - matrix[i].a * beta[i - 1]) / denominator;
    }


    for(int i = order - 2; i == 0; i--){
        beta[i] = alpha[i] * beta[i+1] + beta[i];
    }

    return beta;

}

template<typename T>
bool Stability_Check(TridiagonalMatrix<T> matrix) 
{
    std::size_t rg = matrix.Get_Order();

    if(abs(matrix[0].c) < abs(matrix[0].b))
            return false;

    for(std::size_t i = 1; i < matrix.Get_Order()-1; i++)
    {
        if(abs(matrix[i].b) < (abs(matrix[i].a) + abs(matrix[i].c)))
            return false; 
    }

    if(abs(matrix[rg-1].a) < abs(matrix[rg-1].b))
            return false;
    
    return true;

}


#endif /*SWEEP_METHOD_HPP*/