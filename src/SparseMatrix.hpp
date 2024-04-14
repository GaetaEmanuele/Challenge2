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
/*template <typename T>
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
*/

    template <typename T>
    class SparseMatrixBase {
    public:
        virtual ~SparseMatrixBase() = default;

        virtual T operator()(std::size_t row, std::size_t col) const = 0;
        virtual T& operator()(std::size_t row, std::size_t col) = 0;
        virtual void compress() const = 0;
        virtual void uncompress() = 0;
        virtual std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const = 0;
        virtual void print() const = 0;
    };
    enum class StorageOrder { RowMajor, ColumnMajor };
    template <typename T, StorageOrder Order>
    class Matrix;

    template <typename T>
    class Matrix<T, StorageOrder::RowMajor> : public SparseMatrixBase<T> {
    private:
        std::map<std::array<std::size_t, 2>, T> elements;
        std::vector<std::vector<T>> compressedMatrix;
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;

    public:
        // Implementazione dei metodi per l'ordinamento RowMajor
        Matrix();
        Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);
        T operator()(std::size_t row, std::size_t col) const override;
        T& operator()(std::size_t row, std::size_t col) override;
        void compress() const override;
        void uncompress() override;
        std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        void print() const override;
    };
    //this struct is used for overload the operator < of the map 
    //that stores the no zero element of the matrix. And so define a ColumMajor sparse matrix
    struct ColumnMajorComparator {
    bool operator()(const std::array<std::size_t, 2>& a, const std::array<std::size_t, 2>& b) const {
        // Ordina per colonna (secondo elemento dell'array)
            return a[1] < b[1] || (a[1] == b[1] && a[0] < b[0]); // Secondo elemento (colonna), poi primo elemento (riga)
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
        std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        void print() const override;
    };  */  


}  // namespace algebra


#include "SparseMatrix_impl.hpp"  // Include the implementation file

#endif // SPARSEMATRIX_HPP