#include "SparseMatrix.hpp"

namespace algebra {

template<ScalarOrComplex T,StorageOrder Order>
void read(Matrix<T, Order>& matrix ,const std::string& file_name){
    std::ifstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Error, Impposible open file: " << file_name << std::endl;
    }

    std::string line;
    std::getline(file, line); // read the first line 
    
    // Check on the format of the first line
    if (line.substr(0, 14) != "%%MatrixMarket") {
        std::cerr << "Eror the file is not in format Matrix Market" << std::endl;
    }

    while (std::getline(file, line) && line[0] == '%') {
        // ignore that line since are unusless
    }
    // Read numbers of rows,columns and non zero elements
    
    std::size_t numNonZero;
    std::string numRowsStr, numColsStr, numNonZeroStr;
    std::istringstream iss(line);
    if(!(iss >> numRowsStr >> numColsStr >> numNonZeroStr)) {
        throw std::runtime_error("Error during the reading");
    }
    auto numRows = std::stoul(numRowsStr);
    auto numCols = std::stoul(numColsStr);
    numNonZero = std::stoul(numNonZeroStr);
    std::cout << numRows<< " "<< numCols<< " "<< numNonZero<<std::endl;
    matrix.resize(numRows,numCols);
    // Read the non zero element
    while(std::getline(file,line)){
        std::istringstream elementStream(line);
        std::string nrow,ncol,val;
        elementStream >> nrow >> ncol >> val;
        std::size_t row,col;
        double value;
        row = std::stoul(nrow);
        col = std::stoul(ncol);
        value = std::stod(val);
        matrix.elements[{row-1,col-1}]= static_cast<T>(value);
    }
    file.close();
}

template<ScalarOrComplex T, StorageOrder Order>
std::vector<T> operator*(const Matrix<T,Order>& matrix, const std::vector<T>& vec){
    if constexpr(Order == algebra::StorageOrder::RowMajor){
        // Create a vector to store the result of matrix-vector multiplication
        std::vector<T> result(matrix.numRows,static_cast<T>(0));
        //Check for the correctness of dimension
        if(vec.size() != matrix.numRows){
            std::cerr<<"Error dimension not coeirent"<<std::endl;
            return result;
        }
        if (matrix.isCompressed) {
            // Iterate over the rows of the compressed matrix
            for (std::size_t i = 0; i < matrix.numRows; ++i) {
                T rowSum = static_cast<T>(0);
    
                // Compute the dot product between the current row of the matrix and the vector 'vec'
                for (std::size_t j=0; j< matrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
                        std:: size_t jj = matrix.compressedMatrix[i][j].first;
                        // Multiply the matrix element value with the corresponding vector element and accumulate
                        rowSum += matrix.compressedMatrix[i][j].second * vec[jj];
                        }
                    // Assign the computed dot product to the corresponding index in the result vector
                    result[i] = rowSum;
                }

            return result;
        }else{
            for(std::size_t i=0;i<matrix.numRows;++i){
                //initialization of the value of the ith componets of the result
                T rowSum = static_cast<T>(0);
                //Use Key to exctract the ith row
                std::array<std::size_t,2> Key = {i,0};
                auto it = matrix.elements.lower_bound(Key);
                Key = {i+1,0};
                auto it_end = matrix.elements.lower_bound(Key);
                //Compute Rowsum
                for(;it != it_end;++it){
                    std::size_t j = it->first[1];
                    rowSum += it->second*vec[j];
                }
                //assignment to the return value
                result[i] = rowSum;
            }
            return result;
        }
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
        // Create a vector to store the result of matrix-vector multiplication
    std::vector<T> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.size() != matrix.numRows){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
    if(!matrix.isCompressed){
        for(std::size_t j=0;j<matrix.numCols;++j){
            auto it = matrix.elements.lower_bound({0,j});
            auto it_end = matrix.elements.lower_bound({0,j+1});
            for(;it!=it_end;++it){
                std::size_t jj = it->first[0];
                result[jj] += it->second* vec[jj];
            }
        }
        return result;
    }else{
        for(std::size_t j=0;j<matrix.numCols;++j){
            auto it = matrix.compressedMatrix[j].cbegin();
            auto it_end = matrix.compressedMatrix[j].cend();
            for(;it != it_end;++it){
                std::size_t i = it->first;
                result[i] += it->second*vec[i];
            }
        }
        return result;
    }
    }
}

template<ScalarOrComplex T, StorageOrder Order>
std::vector<T> operator*(const Matrix<T,Order>& matrix,const Matrix<T,Order>& vec){
    if constexpr (Order == algebra::StorageOrder::RowMajor){
        std::vector<T> result(matrix.numRows,static_cast<T>(0));
        if(vec.numRows != matrix.numRows && vec.numCols >1){
            std::cerr<<"Dimension not valid"<<std::endl;
            return result;
        }
        //for semplicity the vector will be used in the compact format
        auto& mutableVec = const_cast<Matrix<T, StorageOrder::RowMajor>&>(vec);
        if(mutableVec.isCompressed){
            mutableVec.uncompress();
        }
        if(matrix.isCompressed){
            // Iterate over the rows of the compressed matrix
            for (std::size_t i = 0; i < matrix.numRows; ++i) {
                T rowSum = static_cast<T>(0);
                // Compute the dot product between the current row of the matrix and the vector 'vec'
                for (std::size_t j=0; j< matrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
                        std:: size_t jj = matrix.compressedMatrix[i][j].first;
                        // Multiply the matrix element value with the corresponding vector element and accumulate
                        if(mutableVec.elements.find({jj,0})!= mutableVec.elements.end()){
                            rowSum += matrix.compressedMatrix[i][j].second * mutableVec.elements.at({jj,0});
                            }
                        }
                    // Assign the computed dot product to the corresponding index in the result vector
                    result[i] = rowSum;
                }
            return result;
        }else{
            for(std::size_t i=0;i<matrix.numRows;++i){
                //initialization of the value of the ith componets of the result
                T rowSum = static_cast<T>(0);
                //Use Key to exctract the ith row
                std::array<std::size_t,2> Key = {i,0};
                auto it = matrix.elements.lower_bound(Key);
                Key = {i+1,0};
                auto it_end = matrix.elements.lower_bound(Key);
                //Compute Rowsum
                for(;it != it_end;++it){
                    std::size_t j = it->first[1];
                    if(mutableVec.elements.find({j,0}) != mutableVec.elements.end()){
                        rowSum += it->second * mutableVec.elements.at({j,0});
                    }
                }
                //assignment to the return value
                result[i] = rowSum;
            }
            return result;
        }
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
             std::vector<T> result(matrix.numRows, static_cast<T>(0));
    //Check for the correctness of dimension
    if(vec.numRows != matrix.numRows and vec.numCols >1 ){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
    //the vector for semplicity must be stored in RowMajor since it is 
    //a column vector
    auto& mutableVec = const_cast<Matrix<T, StorageOrder::RowMajor>&>(vec);
    if(mutableVec.isCompressed){
        mutableVec.uncompress();
    }
       for(std::size_t j=0;j<matrix.numCols;++j){
            auto it = matrix.elements.lower_bound({0,j});
            auto it_end = matrix.elements.lower_bound({0,j+1});
            for(;it!=it_end;++it){
                std::size_t jj = it->first[0];
                if(mutableVec.elements.find({jj,0}) != mutableVec.elements.end()){
                    result[jj] += it->second * mutableVec.at({jj,j});
                }
            }
        }
        return result;
    }
}

}//name space algebra