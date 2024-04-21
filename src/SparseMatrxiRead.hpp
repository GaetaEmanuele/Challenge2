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

}//name space algebra