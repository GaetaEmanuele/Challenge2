#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <cstddef>      // For std::size_t
#include <initializer_list>  // For std::initializer_list

namespace algebra {

    // Matrix class template for a sparse matrix using std::map for dynamic construction
template <typename T>
    class SparseMatrix {
    private:
        std::map<std::array<std::size_t, 2>, T> elements; //data structure used for storing the matrix in uncompressed way
        std::vector<std::vector<T>> compressedMatrix; //data structore for storing it in the comprossure format
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;

    public:
        SparseMatrix();
        SparseMatrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);

        //void operator()(std::size_t row, std::size_t col, T value);
        T operator()(std::size_t row, std::size_t col) const;
        T& operator()(std::size_t row, std::size_t col);

        void compress() const;
        void uncompress();

        std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const;
        void print() const;
    };


}  // namespace algebra


#include "SparseMatrix_impl.hpp"  // Include the implementation file

#endif // SPARSEMATRIX_HPP