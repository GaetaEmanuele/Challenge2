#include "SparseMatrix.hpp"

namespace algebra{

    template<ScalarOrComplex T>
    Matrix<T,StorageOrder:: ColumnMajor>:: Matrix(std::size_t nrow,std::size_t ncol): numRows(nrow),numCols(ncol),isCompressed(false){}

    // Const operator() to access elements in a compressed or uncompressed matrix
    template <ScalarOrComplex T>
    T Matrix<T,StorageOrder::ColumnMajor>::operator()(std::size_t row, std::size_t col) const {
        if (row < numRows && col < numCols) {
                if (isCompressed) {
                    // Access element in compressed matrix
                    auto it = std::find_if(compressedMatrix[col].cbegin(), compressedMatrix[col].cend(),
                            [row](const std::pair<std::size_t, T>& p) {
                                return p.first == row;
                            });
                    if(it != compressedMatrix[col].cend()){
                        return it->second;
                    }
                    throw std::out_of_range("Element not found");
                } else {
                    // Access element in uncompressed map
                    auto it = elements.find({row, col});
                    if (it != elements.end()) {
                        return it->second;
                    } else {
                        throw std::out_of_range("element not found");
                    }
                }
            } else {
                throw std::out_of_range("Index out of boundary");
            }
    }
     // Non-const operator() to modify elements in a compressed or uncompressed matrix
    template <ScalarOrComplex T>
    T& Matrix<T,StorageOrder::ColumnMajor>::operator()(std::size_t row, std::size_t col) {
        if (isCompressed) {
        // Check if the element exists and is non-zero in the compressed matrix
            auto it = std::find_if(compressedMatrix[col].begin(), compressedMatrix[col].end(),
                            [row](const std::pair<std::size_t, T>& p) {
                                return p.first == row;
                            });
            if (row < numRows && col < numCols && it->second != 0) {
                return it->second;
            } else {
                throw std::out_of_range("Indexes out of boundary");
            }
        } else {
             // Access or insert element in uncompressed map
            return elements[{row, col}];
        }
    }
    // Function to compress the sparse matrix representation (const-correct version)
    template <ScalarOrComplex T>
    void Matrix<T,StorageOrder::ColumnMajor>::compress(){
        //The change of form must be done only if the Matrix is not already compress
        if(!isCompressed){
            //Resize the matrix
            compressedMatrix.resize(numCols);
            for(std::size_t i=0;i<numCols;++i){
                auto cit = elements.lower_bound({0,i});
                auto cend = elements.lower_bound({0,i+1});
                for (; cit != cend; ++cit){
                    compressedMatrix[i].emplace_back(cit->first[0],cit->second);
                }
            }
       //Clear the map to free the memory
       elements.clear();

        //Set the compressed flag
        isCompressed = true;
        }
    }

    // Function to uncompress the sparse matrix representation
    template <ScalarOrComplex T>
    void Matrix<T,StorageOrder::ColumnMajor>::uncompress() {
        if (!isCompressed) {
                std::cerr << "Error: Matrix is not compressed. Cannot uncompress." << std::endl;
                return;
            }

            // Rebuild the elements map from compressed matrix
            for(std::size_t i=0;i<numCols;++i){
                for(auto it = compressedMatrix[i].begin();it!=compressedMatrix[i].end();++it){
                    auto jj = it->first;
                    elements[{jj,i}] = it->second;
                }
            }

            // Clear the compressed matrix to free memory
            compressedMatrix.clear();

            // Reset the compression flag
            isCompressed = false;
    }
    // Function to print the matrix (supports both compressed and uncompressed)
    template <ScalarOrComplex T>
    void Matrix<T,StorageOrder::ColumnMajor>::print() const {
        if(isCompressed){
            for(std::size_t i=0;i<numRows;++i){
                for(std::size_t j=0;j<numCols;++j){
                    auto it = std::find_if(compressedMatrix[j].cbegin(), compressedMatrix[j].cend(),
                            [i](const std::pair<std::size_t, T>& p) {
                                return p.first == i;
                            });
                    if(it!=compressedMatrix[i].cend()){
                        std::cout<<it->second;
                    }else{std::cout<<0;}
                }
                std::cout<<std::endl;
            }
        }else{
            for(std::size_t i=0;i<numRows;++i){
                for(std::size_t j=0;j<numCols;++j){
                    std::array<std::size_t,2> Key = {i,j};
                    auto it = elements.find(Key);
                    if(it != elements.cend()){
                        std::cout<<it->second;
                    }else{std::cout<<0;}
                }
                std::cout<<std::endl;
            }
        }
    }
    
    template<ScalarOrComplex T>
    T Matrix<T,StorageOrder::ColumnMajor>::norm(const algebra:: Typenorm& norm_)const{
        if(norm_ == algebra::Typenorm::One){
                //initialization of sum
                T sum = static_cast<T>(0);
                if(isCompressed){
                    //accumulate method + lambda function allows to evaluate the sum for column
                    sum = std::accumulate(compressedMatrix[0].cbegin(),compressedMatrix[0].cend(),static_cast<T>(0),
                                         [](const T& acc,const std::pair<std::size_t,T>&entry){
                                            return acc + std::abs(entry.second);
                                         });
                    for(std::size_t i=1;i<compressedMatrix.size();++i){
                        sum = std::max(sum,std::accumulate(compressedMatrix[i].cbegin(),compressedMatrix[i].cend(),static_cast<T>(0),
                                         [](const T& acc,const std::pair<std::size_t,T>&entry){
                                            return acc + std::abs(entry.second);
                                         }));
                    return sum;
                    }
                }else{
                    //since the Matrix is stored with Column Major order
                    //i exctract the ith column with the map's method lower_bound
                    std::array<std::size_t,2> Key = {0,0};
                    auto it = elements.lower_bound(Key);
                    Key = {0,1};
                    auto it_end = elements.lower_bound(Key);
                    sum = std::accumulate(it,it_end,static_cast<T>(0),
                                         [](const T& acc, const std::pair<const std::array<std::size_t,2>,T>& entry){
                                            return acc + std::abs(entry.second);
                                         }); 
                    for (std::size_t i=1;i<numCols;++i){
                        Key = {0,i};
                        it = elements.lower_bound(Key);
                        Key = {0,i+1};
                        it_end = elements.lower_bound(Key);
                        sum  = std:: max(sum,std::accumulate(it,it_end,static_cast<T>(0),
                                         [](const T& acc, const std::pair<const std::array<std::size_t,2>,T>& entry){
                                            return acc + std::abs(entry.second);
                                         }));
                    }
                    return sum;
                }   
        }else if(norm_ == algebra::Typenorm::Infinity){
                if(isCompressed){
                //StorageOrder = ColumnMajor so i need an auxiliary vector to store
                //the sum by row
                std::vector<T> RowSum(numRows,static_cast<T>(0));
                for(std::size_t i=0;i<numCols;++i){
                    for(std::size_t j=0;j<numRows;++j){
                         //finds the lement (i,j)
                         auto it = std::find_if(compressedMatrix[i].begin(), compressedMatrix[i].end(),
                            [j](const std::pair<std::size_t, T>& p) {
                                return p.first == j;
                            });
                         if (it!= compressedMatrix[i].cend()){
                            RowSum[j] += std::abs(it->second);
                         }
                    }
                }
                return *std::max_element(RowSum.cbegin(),RowSum.cend());
                }else{
                //initialization of the value;
                T sum(static_cast<T>(0)),value;
                value = sum;
                for(std::size_t i=0;i<numRows;++i){
                    //sum if the condition described by the lambda function is true
                    sum = std::accumulate(elements.cbegin(),elements.cend(),static_cast<T>(0),
                            [i](const T& acc, const std::pair<const std::array<size_t, 2>, T>& entry){
                                auto index = entry.first;
                                if(index[0]==i){
                                    return acc + std::abs(entry.second);
                                }
                                return acc;
                            });
                    value = std::max(value,sum);
                } 
                return value;
                }
        }else if(norm_ == algebra::Typenorm::Frobenius){
                T sum(static_cast<T>(0));
                if(isCompressed){
                    T sum = std::accumulate(compressedMatrix[0].cbegin(),compressedMatrix[0].cend(),static_cast<T>(0),[](const T& acc,const std::pair<std::size_t, T>& entry){
                        return acc + std::abs(entry.second*entry.second);});
                    for(std::size_t i=1;i<compressedMatrix.size();++i){
                        sum += std::accumulate(compressedMatrix[i].cbegin(),compressedMatrix[i].cend(),static_cast<T>(0),[](const T& acc,const std::pair<std::size_t, T>& entry){
                            return acc + std::abs(entry.second*entry.second);});
               }
               return sum;
                }else{
                std::array<std::size_t,2> Key = {0,0};
                auto it = elements.lower_bound(Key);
                Key = {0,1};
                auto it_end = elements.lower_bound(Key);
                T sum =std::accumulate(it,it_end,static_cast<T>(0),[](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry){
                    return acc + std::abs(entry.second*entry.second);});
                for(std::size_t i=1;i<numRows;++i){
                    Key = {0,i};
                    it = elements.lower_bound(Key);
                    Key = {0,i+1};
                    it_end = elements.lower_bound(Key);
                    sum +=std::accumulate(it,it_end,static_cast<T>(0),
                                      [](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry){return acc + std::abs(entry.second*entry.second);}
                                      );
                    }
                return sum;
        }
            }
            } 

    template<ScalarOrComplex T>
    void Matrix<T,StorageOrder::ColumnMajor>::resize(std::size_t nrow,std::size_t ncol){
            numRows =nrow;
            numCols = ncol;
            isCompressed =false;
    }
   
} //name space algebra