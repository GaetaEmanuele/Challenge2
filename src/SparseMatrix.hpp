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
    class SparseMatrixBase {
    public:
        virtual ~SparseMatrixBase() = default;

        virtual T operator()(std::size_t row, std::size_t col) const = 0;
        virtual T& operator()(std::size_t row, std::size_t col) = 0;
        virtual void compress() const = 0;
        virtual void uncompress() = 0;
        virtual void print() const = 0;
    };
    
    enum class StorageOrder { RowMajor, ColumnMajor };
    template <typename T, StorageOrder Order>
    class Matrix;

    
    template <typename T>
    std::vector<T> operator*(const Matrix<T, StorageOrder::RowMajor>& matrix, const std::vector<T>& vec);
    template <typename T>
    class Matrix<T, StorageOrder::RowMajor> : public SparseMatrixBase<T> {
    private:
        std::map<std::array<std::size_t, 2>, T> elements;
        std::vector<std::vector<T>> compressedMatrix;
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;

    public:
        // Declaration of the class specification RowMajor
        Matrix();
        Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);
        T operator()(std::size_t row, std::size_t col) const override;
        T& operator()(std::size_t row, std::size_t col) override;
        void compress() const override;
        void uncompress() override;
        friend std::vector<T> algebra::operator*<>(const Matrix<T, StorageOrder::RowMajor>& matrix, const std::vector<T>& vec);
        //std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        void print() const override;
    };
    //this struct is used for overload the operator < of the map 
    //that stores the no zero element of the matrix. And so define a ColumMajor sparse matrix
    struct ColumnMajorComparator {
    bool operator()(const std::array<std::size_t, 2>& a, const std::array<std::size_t, 2>& b) const {
        // now i want to ordinate by column
            return a[1] < b[1] || (a[1] == b[1] && a[0] < b[0]); 
        }
    };

/*

    template <typename T>
    class Matrix<T, StorageOrder::ColumnMajor> : public SparseMatrixBase<T> {
    private:
        std::map<std::array<std::size_t, 2>, T, ColumnMajorComparator> elements;
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;

    public:
        // Implementazione dei metodi per l'ordinamento ColumnMajor
        Matrix();
        Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);
        T operator()(std::size_t row, std::size_t col) const override;
        T& operator()(std::size_t row, std::size_t col) override;
        void compress() const override;
        void uncompress() override;
        friend std::vector<T> operator*(const Matrix<T, StorageOrder::ColumnMajor>& matrix, const std::vector<T>& vec);
        //std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        void print() const override;
    };  */


}  // namespace algebra


#include "SparseMatrix_impl.hpp"  // Include the implementation file
#endif // SPARSEMATRIX_HPP