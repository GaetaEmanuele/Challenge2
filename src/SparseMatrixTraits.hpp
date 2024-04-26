#ifndef SparseMatrixTraits_HPP
#define SparseMatrixTraits_HPP
#include<complex>
#include <iostream>
#include <type_traits>

namespace algebra{
//Concept for undesrtand if a datum is numeric
template<typename T>
concept Numeric = std::is_arithmetic_v<T>;
//Concept for understand if a datum is complex
template<typename T>
concept Complex = requires(T t) {
    { t.real() } -> Numeric;
    { t.imag() } -> Numeric;
};

//Now the type of the elements of the matrix will be or a real/integer value or also 
//a complex value
template<typename T>
concept ScalarOrComplex = Numeric<T> || Complex<T>;
}
#endif /* SparseMatrixTraits_HPP */