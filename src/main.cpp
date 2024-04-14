#include <iostream>
#include "SparseMatrix.hpp"  // Include the header file with SparseMatrix implementation

int main() {
    // Create a sparse matrix with initial non-zero elements
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> mat{
        {0, 0, 1.2},
        {1, 2, -3.4},
        {2, 1, 5.6}
    };

    // Print the initial matrix (uncompressed)
    std::cout << "Initial Matrix (Uncompressed):\n";
    mat.print();
    std::cout << std::endl;

    // Access and modify elements using non-const operator() (uncompressed)
    mat(1, 2) = 10.0;   // Modify existing non-zero element
    mat(2, 0) = 7.8;    // Add a new element (uncompressed)

    // Print the modified matrix (uncompressed)
    std::cout << "Modified Matrix (Uncompressed):\n";
    mat.print();
    std::cout << std::endl;

    // Compress the matrix
    mat.compress();

    // Access and modify elements using non-const operator() (compressed)
    mat(2, 1) = -2.5;   // Modify existing non-zero element
    // Trying to add a new element in compressed state will give an error
    // mat(0, 1) = 3.14;  // Error: Cannot add new elements in compressed matrix

    // Print the modified matrix (compressed)
    std::cout << "Modified Matrix (Compressed):\n";
    mat.print();
    std::cout << std::endl;

    // Uncompress the matrix to add new elements
    mat.uncompress();

    // Add new elements after uncompressing
    mat(0, 1) = 3.14;  // Add a new element after uncompressing

    // Print the final matrix (uncompressed)
    std::cout << "Final Matrix (Uncompressed):\n";
    mat.print();
    std::cout << std::endl;

    // Perform matrix-vector multiplication
    std::vector<double> vec{1.0, 2.0, 3.0};
    auto result = mat.matrixVectorProduct(vec);

    // Print the result of matrix-vector product
    std::cout << "Result of Matrix-Vector Product:\n";
    for (const auto& val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}
