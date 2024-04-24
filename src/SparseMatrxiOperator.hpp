#include "SparseMatrix.hpp"

namespace algebra {

template<typename T,StorageOrder Order>
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

template<typename T, StorageOrder Order>
std::vector<T> operator*(const Matrix<T,Order>& matrix, const std::vector<T>& vec){
    if constexpr(Order == algebra::StorageOrder::RowMajor){
            //The only way to change the state of a const object 
        auto& mutableMatrix = const_cast<Matrix<T, StorageOrder::RowMajor>&>(matrix);
        // Create a vector to store the result of matrix-vector multiplication
        std::vector<T> result(matrix.numRows,static_cast<T>(0));
        //Check for the correctness of dimension
        if(vec.size() != matrix.numRows){
            std::cerr<<"Error dimension not coeirent"<<std::endl;
            return result;
        }
        if (!matrix.isCompressed) {
            mutableMatrix.compress();
        }
        // Iterate over the rows of the compressed matrix
        for (std::size_t i = 0; i < matrix.numRows; ++i) {
            T rowSum = static_cast<T>(0);

            // Compute the dot product between the current row of the matrix and the vector 'vec'
            for (std::size_t j=0; j< mutableMatrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
                std:: size_t jj = mutableMatrix.compressedMatrix[i][j].first;
                // Multiply the matrix element value with the corresponding vector element and accumulate
                rowSum += mutableMatrix.compressedMatrix[i][j].second * vec[jj];
         }

        // Assign the computed dot product to the corresponding index in the result vector
        result[i] = rowSum;
    }

    return result;
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
        // Create a vector to store the result of matrix-vector multiplication
    std::vector<T> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.size() != matrix.numRows){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }

    for(std::size_t j=0;j<matrix.numCols;++j){
        auto it = matrix.elements.lower_bound({0,j});
        auto it_end = matrix.elements.lower_bound({0,j+1});
        for(;it!=it_end;++it){
            std::size_t jj = it->first[0];
            result[jj] += it->second* vec[jj];
        }
    }

    return result;
    }
}

template<typename T, StorageOrder Order>
std::vector<T> operator*(const Matrix<T,Order>& matrix,const Matrix<T,Order>& vec){
    if constexpr (Order == algebra::StorageOrder::RowMajor){
        //The only way to change the state of a const object 
        auto& mutableMatrix = const_cast<Matrix<T, StorageOrder::RowMajor>&>(matrix); 
        auto& mutableVec = const_cast<Matrix<T, StorageOrder::RowMajor>&>(vec);
        // Create a vector to store the result of matrix-vector multiplication
        std::vector<T> result(matrix.numRows, 0);
        //Check for the correctness of dimension
        if(vec.numRows != matrix.numRows and vec.numCols!=1){
            std::cerr<<"Error dimension not coeirent"<<std::endl;
            return result;
        }
        if (!matrix.isCompressed) {
            mutableMatrix.compress();
        }
        if(vec.isCompressed){
            mutableVec.uncompresse();
        }
        // Iterate over the rows of the compressed matrix
        for (std::size_t i = 0; i < matrix.numRows; ++i) {
            T rowSum = static_cast<T>(0);

            // Compute the dot product between the current row of the matrix and the vector 'vec'
            for (std::size_t j=0; j< mutableMatrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
                std:: size_t jj = mutableMatrix.compressedMatrix[i][j].first;
                // Multiply the matrix element value with the corresponding vector element and accumulate
                //i'm using the vec that now is a type of sparse matrix in uncompressed format
                //Since it is simpler find the corrisponding value and the searching operation cost less in a map than in a vector
                if(mutableVec.elements.find({j,jj})!= mutableVec.elements.end()){
                rowSum += matrix.compressedMatrix[i][j].second * mutableVec.elements.at({j,jj});
                }
            }

        // Assign the computed dot product to the corresponding index in the result vector
        result[i] = rowSum;
        }

        return result;
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
             std::vector<T> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.numRows() != matrix.numRows and vec.numCols() !=1 ){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
       for(std::size_t j=0;j<matrix.numCols;++j){
            auto it = matrix.elements.lower_bound({0,j});
            auto it_end = matrix.elements.lower_bound({0,j+1});
            for(;it!=it_end;++it){
                std::size_t jj = it->first[0];
                if(vec.elements.find({jj,j}) != vec.elements.end()){
                    result[jj] += it->second * vec.at({jj,j});
                }
            }
        }
        return result;
    }
}

template<typename T, StorageOrder Order>
std::vector<std::complex<T>> operator*(const Matrix<std::complex<T>,Order>& matrix,const std::vector<std::complex<T>>& vec ){
    if constexpr (Order == algebra::StorageOrder::RowMajor){
             //The only way to change the state of a const object 
    auto& mutableMatrix = const_cast<Matrix<std::complex<T>, StorageOrder::RowMajor>&>(matrix);
   // Create a vector to store the result of matrix-vector multiplication
    std::vector<std::complex<T>> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.size() != matrix.numRows){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
    if (!matrix.isCompressed) {
        mutableMatrix.compress();
    }
    // Iterate over the rows of the compressed matrix
    for (std::size_t i = 0; i < matrix.numRows; ++i) {
        std::complex<T> rowSum = std::complex<T>(0.0, 0.0);

        // Compute the dot product between the current row of the matrix and the vector 'vec'
        for (std::size_t j=0; j< mutableMatrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
            std:: size_t jj = mutableMatrix.compressedMatrix[i][j].first;
            // Multiply the matrix element value with the corresponding vector element and accumulate
            rowSum.real += mutableMatrix.compressedMatrix[i][j].second.real * vec[jj].real - mutableMatrix.compressedMatrix[i][j].second.imag*vec[jj].imag;
            rowSum.imag += mutableMatrix.compressedMatrix[i][j].second.real * vec[jj].imag + mutableMatrix.compressedMatrix[i][j].second.imag * vec[jj].real;
        }

        // Assign the computed dot product to the corresponding index in the result vector
        result[i] = rowSum;
    }

    return result;
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
             // Compress the matrix if it's not already compressed
       //The only way to change the state of a const object 
   // Create a vector to store the result of matrix-vector multiplication
    std::vector<std::complex<T>> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.size() != matrix.numRows){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }

    for(std::size_t j=0;j<matrix.numCols;++j){
        auto it = matrix.elements.lower_bound({0,j});
        auto it_end = matrix.elements.lower_bound({0,j+1});
        for(;it!=it_end;++it){
            std::size_t jj = it->first[0];
            result[jj].real += it->second.real* vec[jj].real - it->second.imag * vec[jj].imag;
            result[jj].imag += it->second.real* vec[jj].imag + it->second.imag * vec[jj].real;      
        }
    }

    return result;
    }
}

template<typename T, StorageOrder Order>
std::vector<std::complex<T>> operator*(const Matrix<std::complex<T>,Order>& matrix, const Matrix<std::complex<T>,Order>& vec){
    if constexpr (Order == algebra::StorageOrder::RowMajor){
             //The only way to change the state of a const object 
    auto& mutableMatrix = const_cast<Matrix<std::complex<T>, StorageOrder::RowMajor>&>(matrix); 
    auto& mutableVec = const_cast<Matrix<std::complex<T>, StorageOrder::RowMajor>&>(vec);
    // Create a vector to store the result of matrix-vector multiplication
    std::vector<std::complex<T>> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.numRows != matrix.numRows and vec.numCols!=1){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
    if (!matrix.isCompressed) {
        mutableMatrix.compress();
    }
    if(vec.isCompressed){
        mutableVec.uncompresse();
    }
    // Iterate over the rows of the compressed matrix
    for (std::size_t i = 0; i < matrix.numRows; ++i) {
        std::complex<T> rowSum = std::complex<T> (0.0,0.0);

        // Compute the dot product between the current row of the matrix and the vector 'vec'
        for (std::size_t j=0; j< mutableMatrix.compressedMatrix[i].size(); ++j) {         // Value of the matrix element
            std:: size_t jj = mutableMatrix.compressedMatrix[i][j].first;
            // Multiply the matrix element value with the corresponding vector element and accumulate
            //i'm using the vec that now is a type of sparse matrix in uncompressed format
            //Since it is simpler find the corrisponding value and the searching operation cost less in a map than in a vector
            if(mutableVec.elements.find({j,jj})!= mutableVec.elements.end()){
            //rowSum += matrix.compressedMatrix[i][j].second * vec.elements.at({j,jj});
              rowSum.real += mutableMatrix.compressedMatrix[i][j].second.real * mutableVec.elements.at({j,jj}).real - mutableMatrix.compressedMatrix[i][j].second.imag*mutableVec.elements.at({j,jj}).imag;
              rowSum.imag += mutableMatrix.compressedMatrix[i][j].second.real * mutableVec.elements.at({j,jj}).imag + mutableMatrix.compressedMatrix[i][j].second.imag * mutableVec.elements.at({j,jj}).real;
            }
        }

        // Assign the computed dot product to the corresponding index in the result vector
        result[i] = rowSum;
    }

    return result;
    }else if constexpr(Order == algebra::StorageOrder::ColumnMajor){
        std::vector<std::complex<T>> result(matrix.numRows, 0);
    //Check for the correctness of dimension
    if(vec.numRows() != matrix.numRows and vec.numCols() !=1 ){
        std::cerr<<"Error dimension not coeirent"<<std::endl;
        return result;
    }
       for(std::size_t j=0;j<matrix.numCols;++j){
            auto it = matrix.elements.lower_bound({0,j});
            auto it_end = matrix.elements.lower_bound({0,j+1});
            for(;it!=it_end;++it){
                std::size_t jj = it->first[0];
                if(vec.elements.find({jj,j}) != vec.elements.end()){
                    result[jj].real += it->secon.real * vec.at({jj,j}).real - it->second.imag*vec.at({jj,j}).imag;
                    result[jj].imag += it->secon.real * vec.at({jj,j}).imag + it->second.imag*vec.at({jj,j}).real;
                }
            }
        }
        return result;
    }
}

}//name space algebra