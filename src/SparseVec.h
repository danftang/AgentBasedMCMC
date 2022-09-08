//
// Created by daniel on 09/04/2021.
//

#ifndef ABMCMC_SPARSEVEC_H
#define ABMCMC_SPARSEVEC_H

#include <vector>
#include <ostream>
#include <algorithm>

template<class T>
class SparseVec {
public:
    std::vector<int> indices;   // array of indices of non-zero elements
    std::vector<T> values; // array of values of non-zero elements

    SparseVec() { }

    explicit SparseVec(int sparseSize): indices(sparseSize), values(sparseSize) { }

    SparseVec(SparseVec<T> &&rvalue) noexcept: indices(std::move(rvalue.indices)), values(std::move(rvalue.values)) { }

    SparseVec(const SparseVec<T> &lvalue): indices(lvalue.indices), values(lvalue.values) { // copy semantics
    }

    explicit SparseVec(const std::vector<T> &dense): indices(), values() {
        for(int i=0; i < dense.size(); ++i) {
            if(double v = dense[i]; v != 0.0) {
                indices.push_back(i);
                values.push_back(v);
            }
        }
    }

    double operator [](int denseIndex) const;

    int sparseSize() const { return indices.size(); }

    std::vector<T> toDense() const;
    std::vector<T> toDense(int dimension) const;
    void insert(int i, T v);
    void clear();
    void resize(size_t size) { indices.resize(size); values.resize(size); }
    void reserve(size_t n) { indices.reserve(n); values.reserve(n);}

    int maxNonZeroIndex() const;


    T operator[](int index) {
        for(int i = 0;i<sparseSize(); ++i) {
            if(indices[i] == index) return values[i];
        }
        return 0.0;
    }

    template<class OTHERNUMBER>
    T operator *(const std::vector<OTHERNUMBER> &other) const {
        T dotProd = 0;
        for(int i=0; i < sparseSize(); ++i) {
            dotProd += values[i] * other[indices[i]];
        }
        return dotProd;
    }

    friend T operator *(const std::vector<T> &lhs, const SparseVec<T> &rhs) {
        return rhs * lhs;
    }

    // copy semantics
    SparseVec<T> &operator =(const SparseVec<T> &lvalue) {
        indices = lvalue.indices;
        values = lvalue.values;
        return *this;
    }

    // move semantics
    SparseVec<T> &operator =(SparseVec<T> &&rvalue) {
        indices = std::move(rvalue.indices);
        values = std::move(rvalue.values);
        return *this;
    }

protected:
    void swap(SparseVec<T> &rvalue) {
        indices.swap(rvalue.indices);
        values.swap(rvalue.values);
    }

};


template<class T>
int SparseVec<T>::maxNonZeroIndex() const {
    int max = 0;
    for(int i=0; i<sparseSize(); ++i) {
        if(indices[i] > max) max = indices[i];
    }
    return max;
}


template<class T>
double SparseVec<T>::operator[](int denseIndex) const {
    auto sparseIterator = std::find(indices.begin(), indices.end(), denseIndex);
    if(sparseIterator != indices.end()) return values[sparseIterator - indices.begin()];
    return 0.0;
}

// Adds element,
// doesn't delete if ind already exists
template<class T>
void SparseVec<T>::insert(int i, T v) {
    if (v != 0.0) {
        indices.push_back(i);
        values.push_back(v);
    }
}


template<class T>
void SparseVec<T>::clear() {
    indices.clear();
    values.clear();
}

// to zero-based dense array
template<class T>
std::vector<T> SparseVec<T>::toDense() const {
    int vecSize = maxNonZeroIndex() + 1;
    std::vector<T> denseVec(vecSize,0);
    for (int i = 0; i < sparseSize(); ++i) {
        denseVec[indices[i]] = values[i];
    }
    return denseVec;
}

template<class T>
std::vector<T> SparseVec<T>::toDense(int dimension) const {
    std::vector<T> denseVec(dimension,0);
    for (int i = 0; i < sparseSize(); ++i) {
        denseVec[indices[i]] = values[i];
    }
    return denseVec;
}


template<class T>
std::ostream &operator<<(std::ostream &out, const SparseVec<T> &sVector) {
    for(auto entry: sVector) {
        out << "[" << entry.index << "] -> " << entry.getValue << "  ";
    }
    return out;
}



#endif //GLPKTEST_SPARSEVEC_H
