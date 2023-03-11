#ifndef CSR_MATRIX
#define CSR_MATRIX

#include <vector>
#include <iostream>
using std::vector;
using std::size_t;


namespace struct_DOC{

    template <typename T>
    struct DOC
    {
        size_t i, j;
        T value;
    };
}

template <typename T>
class CsrMatrix{
private:

    vector<T> elems;
    vector<size_t> col_ind;
    vector<size_t> amount_elems;

    size_t size_str; 
    size_t size_col;

public:
    
    CsrMatrix(const vector<struct_DOC::DOC<T>> &vector_of_matrix)
    {
        elems.resize(vector_of_matrix.size());
        col_ind.resize(vector_of_matrix.size());
        amount_elems.resize(vector_of_matrix[vector_of_matrix.size() - 1].i + 2);

        amount_elems[0] = 0;

        size_t j = 1;

        size_col = 0;

        for (size_t it = 0; it < vector_of_matrix.size(); it++)
        {
            elems[it]   = vector_of_matrix[it].value;
            
            col_ind[it] = vector_of_matrix[it].j;
            
            if (col_ind[it] > size_col){
                size_str = col_ind[it];
            }

            
            if (it > 0 && (vector_of_matrix[it].i - vector_of_matrix[it - 1].i != 0))
            {
                for (size_t k = 0; k < (vector_of_matrix[it].i - vector_of_matrix[it - 1].i); k++)
                {
                    amount_elems[j] = it;
                    j++;
                }
            }
        }

        amount_elems[j] = vector_of_matrix.size();

        size_col = amount_elems.size() - 1;
        size_str++;

    }


    size_t SizeStr() const {
        return size_str;
    }


    size_t SizeCol() const {
        return size_col;
    }


    const T operator() (size_t i, size_t j) const{

    for(size_t k = amount_elems[i]; k < amount_elems[i + 1]; k++)
    {
        if(j == col_ind[k])
            return elems[k];
    }

    return 0;
}


    std::vector<T> operator*(const std::vector<T> &other_vec) const {

        std::vector<T> result_vec;
        result_vec.reserve(size_col);

        for (size_t i = 0; i < this->size_col; i++){
            
            T result_elem = 0;
            
            for (size_t j = this->amount_elems[i]; j < amount_elems[i + 1]; j++) {
                result_elem += other_vec[this->col_ind[j]] * this->elems[j];
                
                
            }

            result_vec.push_back(result_elem);

        }

        return result_vec;
    }


    size_t GetAmountNozeroElems() const {
        return amount_elems[amount_elems.size() - 1];
    }

    CsrMatrix<T> GetDownTriagMatr() const{

        vector<struct_DOC::DOC<T>> down_vec;
        down_vec.reserve(elems.size());

        for(size_t i = 0; i < size_col; i++){
            for(size_t j = amount_elems[i]; j < amount_elems[i + 1]; j++){
                if(i > col_ind[j]){
                    down_vec.push_back({i, col_ind[j], elems[j]});
                }
            }
        }

        CsrMatrix result_matr(down_vec);
        return result_matr;
    }

    CsrMatrix<T> GetUpTriagMatr() const{

        vector<struct_DOC::DOC<T>> up_vec;
        up_vec.reserve(elems.size());

        for(size_t i = 0; i < size_col; i++){
            for(size_t j = amount_elems[i]; j < amount_elems[i + 1]; j++){
                if(i < col_ind[j]){
                    up_vec.push_back({i, col_ind[j], elems[j]});
                }
            }
        }

        CsrMatrix result_matr(up_vec);
        return result_matr;
    }

    CsrMatrix<T> GetDiagMatr() const{

        vector<struct_DOC::DOC<T>> diag_vec;
        diag_vec.reserve(elems.size());

        for(size_t i = 0; i < size_col; i++){
            for(size_t j = amount_elems[i]; j < amount_elems[i + 1]; j++){
                if(i == col_ind[j]){
                    diag_vec.push_back({i, i, elems[j]});
                }
            }
        }

        CsrMatrix result_matr(diag_vec);
        return result_matr;
    }

    CsrMatrix<T> operator+(const CsrMatrix<T>&  other_matr) const {
        
        size_t N = this->GetAmountNozeroElems() + other_matr.GetAmountNozeroElems();
        
        vector<struct_DOC::DOC<T>> vec_of_sum;
        vec_of_sum.reserve(N);

        for(size_t i = 0; i < size_col; i++){
        
            for(size_t j = 0; j < size_str; j++){
                if(((*this)(i,j) + other_matr(i,j)) > 0){
                    vec_of_sum.push_back({i, j, (*this)(i,j) + other_matr(i,j)});
                }

                else continue;
            }
        
        }

        CsrMatrix result_matr(vec_of_sum);
        return result_matr;
    }


};






#endif
