#ifndef SPARSEMATRIX_HPP
#define SPARSEMATRIX_HPP

#include <iostream>
#include <map>
#include <array>
#include <vector>
#include <initializer_list>  // For std::initializer_list
#include <complex>
#include <string>
#include <algorithm>
#include<utility>
#include<cstddef>
#include<stdexcept>
#include<numeric>
#include<cmath>
#include <fstream>
namespace algebra {

    // Matrix class template for a sparse matrix using std::map for dynamic construction
    enum class Typenorm { One, Infinity,Frobenius };
    template <typename T>
    class SparseMatrixBase {
    public:
        virtual ~SparseMatrixBase() = default;

        virtual T operator()(std::size_t row, std::size_t col) const = 0;
        virtual T& operator()(std::size_t row, std::size_t col) = 0;
        virtual T norm(const algebra:: Typenorm&)const=0;
        virtual void compress()= 0;
        virtual void uncompress() = 0;
        virtual void print() const = 0;
    };
    
    enum class StorageOrder { RowMajor, ColumnMajor };
    
    template <typename T, StorageOrder Order>
    class Matrix;

    template <typename T, StorageOrder Order>
    void read(Matrix<T, Order>& matrix, const std::string& file_name);
    
    template <typename T>
    std::vector<T> operator*(const Matrix<T, StorageOrder::RowMajor>& matrix, const std::vector<T>& vec);
    
    template<typename T>
    std::vector<T> operator* (const Matrix<T, StorageOrder::RowMajor>& matrix, const Matrix<T,StorageOrder::RowMajor>& vec);
    
    template <typename T>
    class Matrix<T, StorageOrder::RowMajor> : public SparseMatrixBase<T> {
    private:
        std::map<std::array<std::size_t, 2>, T> elements;
        std::vector<std::vector<std::pair<std::size_t,T>>> compressedMatrix;
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;

    public:
        // Declaration of the class specification RowMajor
        Matrix(std::size_t nrow,std::size_t ncol);
        Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);
        bool is_compressed()const{return isCompressed;};
        T operator()(std::size_t row, std::size_t col) const override;
        void resize(std::size_t nrow,std::size_t ncol);
        T& operator()(std::size_t row, std::size_t col) override;
        void compress() override;
        void uncompress() override;
        T norm(const algebra:: Typenorm& norm_)const override;
        friend std::vector<T> algebra::operator*<>(const Matrix<T, StorageOrder::RowMajor>& matrix, const std::vector<T>& vec);
        friend std::vector<T> algebra::operator*<>(const Matrix<T, StorageOrder::RowMajor>& matrix, const Matrix<T,StorageOrder::RowMajor>& vec);
        //std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        friend void read<>(Matrix<T, StorageOrder::RowMajor>& matrix,const std::string& file_name);
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

    template <typename T>
    std::vector<T> operator*(const Matrix<T, StorageOrder::ColumnMajor>& matrix, const std::vector<T>& vec);
    template <typename T>
    std::vector<T> operator*(const Matrix<T, StorageOrder::ColumnMajor>& matrix, Matrix<T, StorageOrder::ColumnMajor>& vec);

    template<typename T>
    class Matrix<T, StorageOrder::ColumnMajor> : public SparseMatrixBase<T> {
    private:
        std::map<std::array<std::size_t, 2>, T, ColumnMajorComparator> elements;
        std::vector<std::vector<std::pair<std::size_t,T>>> compressedMatrix;
        std::size_t numRows;
        std::size_t numCols;
        bool isCompressed;
    public:
        // Declaration for colum major
        Matrix(std::size_t nrow,std::size_t ncol);
        Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList);
        bool is_compressed()const{return isCompressed;};
        void resize(std::size_t nrow,std::size_t ncol);
        T operator()(std::size_t row, std::size_t col) const override;
        T& operator()(std::size_t row, std::size_t col) override;
        void compress() override;
        void uncompress() override;
        T norm(const algebra::Typenorm& norm_)const override;
        friend void read<>(Matrix<T,StorageOrder::ColumnMajor>& matrix,const std::string& file_name);
        friend std::vector<T> operator*<>(const Matrix<T, StorageOrder::ColumnMajor>& matrix, const std::vector<T>& vec);
        friend std::vector<T> operator*<>(const Matrix<T, StorageOrder::ColumnMajor>& matrix, const Matrix<T,StorageOrder::ColumnMajor>& vec);
        //std::vector<T> matrixVectorProduct(const std::vector<T>& vec) const override;
        void print() const override;
    };  


}  // namespace algebra


#include "SparseMatrix_impl.hpp"  // Include the implementation file for RowMajor
#include "SparseMatrix_impl_col.hpp" //Include the implementation fil for ColumnMajor
#include "SparseMatrxiRead.hpp"
#endif // SPARSEMATRIX_HPP