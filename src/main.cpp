#include <iostream>
#include "SparseMatrix.hpp"  // Include the header file with SparseMatrix implementation
#include "chrono.hpp"
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
    auto result = mat * vec;

    // Print the result of matrix-vector product
    std::cout << "Result of Matrix-Vector Product:\n";
    for (const auto& val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    Timings::Chrono chrono;
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> A(0,0);
    chrono.start();
    algebra::read(A,"lnsp_131.mtx");
    chrono.stop();
    double t = chrono.wallTime();
    std::cout<<"The reading requires: "<<t<<"mics"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"The matrix A is";
    std::cout<<std::endl;
    A.print();
    std::cout<<std::endl;
    std::cout <<"The One norm of A is:";
    double norm1 = A.norm(algebra::Typenorm::One);
    std::cout<<norm1;
    std::cout<<std::endl;
    std::cout <<"The Infinity norm of A is:";
    double norm2 = A.norm(algebra::Typenorm::Infinity);
    std::cout<<norm2;
    std::cout<<std::endl;
    std::cout <<"The Frobenoius norm of A is:";
    double norm3 = A.norm(algebra::Typenorm::Frobenius);
    std::cout<<norm3;
    std::cout<<std::endl;
    std::vector<double> a(131,1.0);
    std::cout<<"The reusult of the matrix*vector multiplication is :"<<std::endl;
    chrono.start();
    std::vector<double> x =  A *a;
    chrono.stop();
    t = chrono.wallTime();
    std::cout<<"The Matrix vector product requires: "<<t<<"micsec"<<"RowMajor"<<std::endl;
    for(std::size_t i=0;i<131;++i){
        std::cout<<x[i];
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> B(0,0);
    algebra::read(B,"lnsp_131.mtx");
    std::cout<<std::endl;
    B.print();
    std::cout<<std::endl;
    std::cout <<"The One norm of B is:";
    norm1 = B.norm(algebra::Typenorm::One);
    std::cout<<norm1;
    std::cout<<std::endl;
    std::cout <<"The Infinity norm of B is:";
    norm2 = B.norm(algebra::Typenorm::Infinity);
    std::cout<<norm2;
    std::cout<<std::endl;
    std::cout <<"The Frobenoius norm of B is:";
    norm3 = B.norm(algebra::Typenorm::Frobenius);
    std::cout<<norm3;
    std::cout<<std::endl;
    std::cout<<"The reusult of the matrix*vector multiplication is :"<<std::endl;
    chrono.start();
    x =  B *a;
    chrono.stop();
    std::cout<<std::endl;
    std::cout<<"the mtrix vector multiplication requires "<<chrono.wallTime()<<" micsec"<<"ColumnMajor"<<std::endl;
    for(std::size_t i=0;i<131;++i){
        std::cout<<x[i];
    }
    return 0; 
}
