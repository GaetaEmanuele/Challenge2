#include <iostream>
#include "SparseMatrix.hpp"  // Include the header file with SparseMatrix implementation
#include "chrono.hpp" //Include of the Professor's Formaggia Utilities
int main() {
    std::cout<<"#test with a sparse matrix stored in RowMajor order#"<<std::endl;
    Timings::Chrono chrono;
    //initialization of the matrix
    algebra::Matrix<double,algebra::StorageOrder::RowMajor> A(0,0);
    chrono.start();
    //reading the matrix from a file
    algebra::read(A,"lnsp_131.mtx");
    chrono.stop();
    double t = chrono.wallTime();
    //print the time spent for reading the matrix
    std::cout<<"The reading requires: "<<t<<" mics"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"The matrix A is";
    std::cout<<std::endl;
    //call of the method print of the class
    A.print();
    std::cout<<std::endl;
    //evaluation of the one norm
    std::cout <<"NORM ONE : ";
    double norm1 = A.norm(algebra::Typenorm::One);
    std::cout<<norm1;
    std::cout<<std::endl;
    //evaluation of the infinity norm
    std::cout <<"NORM INFINITY: ";
    double norm2 = A.norm(algebra::Typenorm::Infinity);
    std::cout<<norm2;
    std::cout<<std::endl;
    //evaluation of the Frobenious norme
    std::cout <<"FROBENIOUS NORME: ";
    double norm3 = A.norm(algebra::Typenorm::Frobenius);
    std::cout<<norm3;
    std::cout<<std::endl;
    //Initialization of the vector
    std::vector<double> a(131,1.0);
    std::cout<<"Matrix, Vector multiplication :"<<std::endl;
    chrono.start();
    std::vector<double> x =  A *a;
    chrono.stop();
    t = chrono.wallTime();
    std::cout<<"The Matrix vector product requires: "<<t<<" micsec"<<std::endl;
    for(std::size_t i=0;i<131;++i){
        std::cout<<x[i];
    }
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"#test with a sparse matrix stored in ColumnMajor order#"<<std::endl;
    //Initialization of the matrix
    algebra::Matrix<double,algebra::StorageOrder::ColumnMajor> B(0,0);
    //Reading the matrix from a file
    algebra::read(B,"lnsp_131.mtx");
    std::cout<<std::endl;
    //call the print method of th class
    B.print();
    std::cout<<std::endl;
    //evaluation of the one norm
    std::cout <<"ONE NORM: ";
    norm1 = B.norm(algebra::Typenorm::One);
    std::cout<<norm1;
    std::cout<<std::endl;
    //evaluation of the infinity norm
    std::cout <<"INFINITY NORM: ";
    norm2 = B.norm(algebra::Typenorm::Infinity);
    std::cout<<norm2;
    std::cout<<std::endl;
    //evaluation of the Frobenious norm
    std::cout <<"FROBENIOUS NORM: ";
    norm3 = B.norm(algebra::Typenorm::Frobenius);
    std::cout<<norm3;
    std::cout<<std::endl;
    std::cout<<"Matrix, Vector multiplication: "<<std::endl;
    chrono.start();
    //matrix vector multiplication
    x =  B *a;
    chrono.stop();
    std::cout<<std::endl;
    std::cout<<"the mtrix vector multiplication requires "<<chrono.wallTime()<<" micsec"<<std::endl;
    //print the result
    for(std::size_t i=0;i<131;++i){
        std::cout<<x[i];
    }
 
    algebra:: Matrix<double,algebra::StorageOrder::RowMajor> vec(131,1);
    for(std::size_t i=0;i<131;++i){
        vec(i,0) = 1.0;
    }
    std::cout<<std::endl;
    std::cout<<"Matrix Column Vector moltupliucation:"<<std::endl;
    std::vector<double> res = A*vec;
    for(std::size_t i=0;i<res.size();++i){
        std::cout<<res[i];
    }
    std::cout<<std::endl;
    
    return 0; 
}
