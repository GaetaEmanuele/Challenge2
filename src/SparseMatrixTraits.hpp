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
}
#endif /* SparseMatrixTraits_HPP */