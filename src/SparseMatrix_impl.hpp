#include "SparseMatrix.hpp"

namespace algebra {
    
    template<ScalarOrComplex T>
    Matrix<T,StorageOrder::RowMajor>:: Matrix(std::size_t nrow,std::size_t ncol): numRows(nrow),numCols(ncol),isCompressed(false){}
    
    template<ScalarOrComplex T>
    void Matrix<T,StorageOrder::RowMajor>::resize(std::size_t nr,std::size_t nc){
        numRows=nr;
        numCols=nc;
        isCompressed = false;
    }
    // Const operator() to access elements in a compressed or uncompressed matrix
    template <ScalarOrComplex T>
    T Matrix<T,StorageOrder::RowMajor>::operator()(std::size_t row, std::size_t col) const {
        if (row < numRows && col < numCols) {
                if (isCompressed) {
                    // Access element in compressed matrix
                        auto it = std::find_if(compressedMatrix[row].cbegin(), compressedMatrix[row].cend(),
                            [col](const std::pair<std::size_t, T>& p) {
                                return p.first == col;
                            });
                        if(it!= compressedMatrix[row].cend()){
                            return it->second;
                        }
                    throw std::out_of_range("Element not exist");
                } else {
                    // Access element in uncompressed map
                    auto it = elements.find({row, col});
                    if (it != elements.end()) {
                        return it->second;
                    } else {
                        return T(); // Element not found, return default value (0 for numeric types)
                    }
                }
            } else {
                throw std::out_of_range("Index out of boundary");
            }
    }
     // Non-const operator() to modify elements in a compressed or uncompressed matrix
    template <ScalarOrComplex T>
    T& Matrix<T,StorageOrder::RowMajor>::operator()(std::size_t row, std::size_t col) {
        if (isCompressed) {
        // Check if the element exists and is non-zero in the compressed matrix
             auto it = std::find_if(compressedMatrix[row].begin(), compressedMatrix[row].end(),
                            [col](const std::pair<std::size_t, T>& p) {
                                return p.first == col;
                            });
            if (row < numRows && col < numCols && it!=compressedMatrix[row].cend()&&it->second != 0) {
                return it->second;
            } else {
                //If the element is not foud this will cause an exception
                throw std::out_of_range("Index out of boundary");
            }
        } else {
             // Access or insert element in uncompressed map
            return elements[{row, col}];
        }
    }
    // Function to compress the sparse matrix representation (const-correct version)
    template <ScalarOrComplex T>
    void Matrix<T,StorageOrder::RowMajor>::compress()  { 
        //The change of form must be done only if the Matrix is not already compress
        if(!isCompressed){
        compressedMatrix.resize(numRows);
        for(std::size_t i=0;i<numRows;++i){
            auto cit = elements.lower_bound({i,0});
            auto cend = elements.lower_bound({i+1,0});
            for (; cit != cend; ++cit){
                compressedMatrix[i].emplace_back(cit->first[1],cit->second);
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
    void Matrix<T,StorageOrder::RowMajor>::uncompress() {
        if (!isCompressed) {
                std::cerr << "Error: Matrix is not compressed. Cannot uncompress." << std::endl;
                return;
            }

            // Rebuild the elements map from compressed matrix
            for(std::size_t i=0;i<numRows;++i){
                for(auto it = compressedMatrix[i].begin();it!=compressedMatrix[i].end();++it){
                    auto jj = it->first;
                    elements[{i,jj}] = it->second;
                }
            }

            // Clear the compressed matrix to free memory
            compressedMatrix.clear();

            // Reset the compression flag
            isCompressed = false;
    }
    // Function to print the matrix (supports both compressed and uncompressed)
    template <ScalarOrComplex T>
    void Matrix<T,StorageOrder::RowMajor>::print() const {
        if(isCompressed){
            for(std::size_t i=0;i<numRows;++i){
                for(std::size_t j=0;j<numCols;++j){
                    auto it = std::find_if(compressedMatrix[i].cbegin(), compressedMatrix[i].cend(),
                            [j](const std::pair<std::size_t, T>& p) {
                                return p.first == j;
                            });
                    if(it!=compressedMatrix[i].cend()){
                        std::cout<<it->second;
                    }else{std::cout<<0;}
                }
                std::cout<<std::endl;
            }
        }else{
            for(std::size_t i=0;i<numRows;++i){
                auto it_b = elements.lower_bound({i,0});
                auto it_e = elements.lower_bound({i+1,0});
                for(std::size_t j=0;j<numCols;++j){
                    auto it = std::find_if(it_b,it_e,
                            [j](const std::pair<std::array<std::size_t,2>,T>& p) {
                                return p.first[1] == j;
                            });
                    if(it != it_e){
                        std::cout<<it->second;
                    }else{std::cout<<0;}
                }
                std::cout<<std::endl;
            }
        }
    }

    template<ScalarOrComplex T>
    T Matrix<T,StorageOrder::RowMajor>::norm(const algebra:: Typenorm& norm_)const{
       if(norm_ == algebra::Typenorm::One){
            if(isCompressed){
                //initialization
                std::vector<T> ColumnSum(numCols,static_cast<T>(0));
                for(std::size_t i=0;i<numRows;++i){
                    for(std::size_t j=0;j<numCols;++j){
                         //find the element(i,j)
                         auto it = std::find_if(compressedMatrix[i].cbegin(), compressedMatrix[i].cend(),
                            [j](const std::pair<std::size_t, T>& p) {
                                return p.first == j;
                            });
                         //do the sum if and only if the value(i,j) was found
                         if (it!= compressedMatrix[i].cend()){
                            ColumnSum[j] += std::abs(it->second);
                         }
                    }
                }
                return *std::max_element(ColumnSum.cbegin(),ColumnSum.cend());
            }else{
                //initialization
                T sum(static_cast<T>(0)),value;
                value = sum;
                for(std::size_t i=0;i<numCols;++i){
                    //do the sum if the condition is satisfies (if the corresponding value was found)
                    sum = std::accumulate(elements.cbegin(),elements.cend(),static_cast<T>(0),
                            [i](const T& acc, const std::pair<const std::array<size_t, 2>, T>& entry){
                                auto index = entry.first;
                                if(index[1]==i){
                                    return acc + std::abs(entry.second);
                                }
                                return acc;
                            });
                    value = std::max(value,sum);
                } 
                return value;
            }
       }else if(norm_ == algebra :: Typenorm::Infinity){
            T sum = static_cast<T>(0);
            if(isCompressed){
               //sum of the first row when the elemet is found
               sum =  std::accumulate(compressedMatrix[0].cbegin(),compressedMatrix[0].cend(),static_cast<T>(0),
               [](const T& acc,const std::pair<std::size_t, T>& entry){
                return acc + std::abs(entry.second);
               });
               for(std::size_t i=1;i<compressedMatrix.size();++i){
                sum = std::max(sum ,std::accumulate(compressedMatrix[i].cbegin(),compressedMatrix[i].cend(),static_cast<T>(0),
               [](const T& acc, const std::pair<std::size_t, T>& entry){
                return acc + std::abs(entry.second);
               }));
               }
               return sum;
            }else{
                std::array<std::size_t,2> Key = {0,0};
                //Order = RowMajor so i can exctract the ith row with low_bound
                auto it = elements.lower_bound(Key);
                Key = {1,0};
                auto it_end = elements.lower_bound(Key);
                sum = std::accumulate(it,it_end,static_cast<T>(0),
                [](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry){
                    return acc + std::abs(entry.second);
                });
                for(std::size_t i=1;i<numRows;++i){
                    Key = {i,0};
                    it = elements.lower_bound(Key);
                    Key = {i+1,0};
                    it_end = elements.lower_bound(Key);
                    sum = std::max(sum,std::accumulate(it,it_end,static_cast<T>(0),
                                      [](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry ){return acc + std::abs(entry.second);}
                                      ));
                }
                return sum;
            }
       }else if(norm_ == algebra::Typenorm::Frobenius){
                if(isCompressed){
                    T sum = std::accumulate(compressedMatrix[0].cbegin(),compressedMatrix[0].cend(),static_cast<T>(0),[](const T& acc,const std::pair<std::size_t, T>& entry){
                        return acc + std::abs(entry.second*entry.second);});
                    for(std::size_t i=1;i<compressedMatrix.size();++i){
                        sum += std::accumulate(compressedMatrix[i].cbegin(),compressedMatrix[i].cend(),static_cast<T>(0),[](const T& acc,const std::pair<std::size_t, T>& entry){
                            return acc + std::abs(entry.second*entry.second);});
               }
               return std::sqrt(sum);
                }else{
                std::array<std::size_t,2> Key = {0,0};
                auto it = elements.lower_bound(Key);
                Key = {1,0};
                auto it_end = elements.lower_bound(Key);
                T sum =std::accumulate(it,it_end,static_cast<T>(0),[](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry){
                    return acc + std::abs(entry.second*entry.second);});
                for(std::size_t i=1;i<numRows;++i){
                    Key = {i,0};
                    it = elements.lower_bound(Key);
                    Key = {i+1,0};
                    it_end = elements.lower_bound(Key);
                    sum +=std::accumulate(it,it_end,static_cast<T>(0),
                                      [](const T& acc,const std::pair<const std::array<size_t, 2>, T>& entry){return acc + std::abs(entry.second*entry.second);}
                                      );
                }
                return std::sqrt(sum);
                }
       }
    }
}  // namespace algebra

