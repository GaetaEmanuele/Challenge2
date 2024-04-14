#include "SparseMatrix.hpp"

namespace algebra {

    template <typename T>
    Matrix<T,StorageOrder::RowMajor>::Matrix() : numRows(0), numCols(0), isCompressed(false) {
        // Constructor without initialize List
    }

    template <typename T>
    Matrix<T,StorageOrder::RowMajor>::Matrix(std::initializer_list<std::tuple<std::size_t, std::size_t, T>> initList)
        : numRows(0), numCols(0), isCompressed(false) {
        // Constructor with initializer list of non zero elem
         for (const auto& tuple : initList) {
                std::size_t row = std::get<0>(tuple);
                std::size_t col = std::get<1>(tuple);
                T value = std::get<2>(tuple);

                // Create a key using the row and column indices
                std::array<std::size_t, 2> key = {row, col};

                // Insert the element into the map
                elements[key] = value;

                // Update the dimensions of the matrix
                numRows = std::max(numRows, row + 1);
                numCols = std::max(numCols, col + 1);
            }
    }
    // Function to insert or update a non-zero element in the matrix
    /*template <typename T>
    void SparseMatrix<T>::operator()(std::size_t row, std::size_t col, T value) {
        // Create a key using the row and column indices
            std::array<std::size_t, 2> key = {row, col};

            // Insert or update the element in the map
            elements[key] = value;

            // Update the dimensions of the matrix
            numRows = std::max(numRows, row + 1);
            numCols = std::max(numCols, col + 1);

            // Reset the compression flag
            isCompressed = false;
    }*/
    // Const operator() to access elements in a compressed or uncompressed matrix
    template <typename T>
    T Matrix<T,StorageOrder::RowMajor>::operator()(std::size_t row, std::size_t col) const {
        if (row < numRows && col < numCols) {
                if (isCompressed) {
                    // Access element in compressed matrix
                    return compressedMatrix[row][col];
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
                std::cerr << "Error: Index out of bounds." << std::endl;
                return T(); // Return default value if index is out of bounds
            }
    }
     // Non-const operator() to modify elements in a compressed or uncompressed matrix
    template <typename T>
    T& Matrix<T,StorageOrder::RowMajor>::operator()(std::size_t row, std::size_t col) {
        if (isCompressed) {
        // Check if the element exists and is non-zero in the compressed matrix
            if (row < numRows && col < numCols && compressedMatrix[row][col] != 0) {
                return compressedMatrix[row][col];
            } else {
                std::cerr << "Error: Cannot add new elements in compressed matrix. Uncompress first to modify." << std::endl;
                // Returning a reference to a temporary, not safe
                // Consider throwing an exception or handling differently based on your needs
                return compressedMatrix[0][0]; // Return a reference to the first element (not safe)
            }
        } else {
             // Access or insert element in uncompressed map
            return elements[{row, col}];
        }
    }
    // Function to compress the sparse matrix representation (const-correct version)
    template <typename T>
    void Matrix<T,StorageOrder::RowMajor>::compress() const {
        if (!isCompressed) {
            auto& self = const_cast<Matrix<T,StorageOrder::RowMajor>&>(*this); // Cast away constness

            // Initialize the compressed matrix with zeros
            self.compressedMatrix.assign(numRows, std::vector<T>(numCols, 0));

            // Copy elements from map to compressed matrix
            for (const auto& elem : self.elements) {
                std::size_t row = elem.first[0];
                std::size_t col = elem.first[1];
                T value = elem.second;
                self.compressedMatrix[row][col] = value;
            }

        // Clear the map to free memory
        self.elements.clear();

        // Set the compression flag
        self.isCompressed = true;
    }
    }

    // Function to uncompress the sparse matrix representation
    template <typename T>
    void Matrix<T,StorageOrder::RowMajor>::uncompress() {
        if (!isCompressed) {
                std::cerr << "Error: Matrix is not compressed. Cannot uncompress." << std::endl;
                return;
            }

            // Rebuild the elements map from compressed matrix
            for (std::size_t i = 0; i < numRows; ++i) {
                for (std::size_t j = 0; j < numCols; ++j) {
                    if (compressedMatrix[i][j] != 0) {
                        elements[{i, j}] = compressedMatrix[i][j];
                    }
                }
            }

            // Clear the compressed matrix to free memory
            compressedMatrix.clear();

            // Reset the compression flag
            isCompressed = false;
    }

    // Function to perform matrix-vector multiplication
    template <typename T>
    std::vector<T> Matrix<T,StorageOrder::RowMajor>::matrixVectorProduct(const std::vector<T>& vec) const {
        if (!isCompressed) {compress();}
                // Perform matrix-vector product using compressed representation
                std::vector<T> result(numRows, 0);
                for (std::size_t i = 0; i < numRows; ++i) {
                    for (std::size_t j = 0; j < numCols; ++j) {
                        result[i] += compressedMatrix[i][j] * vec[j];
                    }
                }
                return result;
    }
    // Function to print the matrix (supports both compressed and uncompressed)
    template <typename T>
    void Matrix<T,StorageOrder::RowMajor>::print() const {
        if (isCompressed) {
                // Print the compressed matrix
                for (std::size_t i = 0; i < numRows; ++i) {
                    for (std::size_t j = 0; j < numCols; ++j) {
                        std::cout << compressedMatrix[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
            } else {
                // Print the elements map (uncompressed)
                for (std::size_t i = 0; i < numRows; ++i) {
                    for (std::size_t j = 0; j < numCols; ++j) {
                        std::array<std::size_t, 2> key = {i, j};
                        auto it = elements.find(key);
                        if (it != elements.end()) {
                            std::cout << it->second << " ";
                        } else {
                            std::cout << "0 ";
                        }
                    }
                    std::cout << std::endl;
                }
    };
    }

}  // namespace algebra
