//
// Created by daniel on 09/04/2021.
//

#ifndef ABMCMC_SPARSEVEC_H
#define ABMCMC_SPARSEVEC_H

#include <vector>
#include <map>
#include <ostream>
#include <algorithm>
#include <boost/serialization/access.hpp>


template<class T>
class SparseVec {
public:
    std::vector<int> indices;   // array of indices of non-zero elements
    std::vector<T> values; // array of values of non-zero elements

    SparseVec()=default;
    SparseVec(const SparseVec<T> &other)=default;
    SparseVec(SparseVec<T> &&other)=default;

    explicit SparseVec(int sparseSize): indices(sparseSize), values(sparseSize) { }

    explicit SparseVec(const std::vector<T> &dense) {
        for(int i=0; i < dense.size(); ++i) {
            if(double v = dense[i]; v != 0.0) {
                indices.push_back(i);
                values.push_back(v);
            }
        }
    }

    explicit SparseVec(const std::map<int,T> &map) {
        reserve(map.size());
        for(const auto &entry: map) insert(entry.first, entry.second);
    }

    SparseVec(const std::initializer_list<std::pair<int,T>> &initValues) {
        reserve(initValues.size());
        for(const std::pair<int,T> &entry: initValues) insert(entry.first, entry.second);
    }

    int sparseSize() const { return indices.size(); }

    std::vector<T> toDense() const;
    std::vector<T> toDense(int dimension) const;
    void insert(int i, const T &v); // Warning: none of the inserts check for duplicate entries
    void insert(int i, T &&v);
    void insert(const SparseVec<T> &other);
    void clear();
    void resize(size_t size) { indices.resize(size); values.resize(size); }
    void reserve(size_t n) { indices.reserve(n); values.reserve(n);}

    int maxNonZeroIndex() const;

    // TODO: need an accessor class for indexing
//    T operator [](int denseIndex) const {
//        auto sparseIterator = std::find(indices.begin(), indices.end(), denseIndex);
//        if(sparseIterator != indices.end()) return values[sparseIterator - indices.begin()];
//        return 0;
//    }
//
//    T operator[](int index) {
//        for(int i = 0;i<sparseSize(); ++i) {
//            if(indices[i] == index) return values[i];
//        }
//        return 0;
//    }

    template<typename ELE, typename = decltype(std::declval<T &>() *= std::declval<const ELE &>())>
    SparseVec &operator *=(const ELE &element) {
        for(int i=0; i < sparseSize(); ++i) values[i] *= element;
        return *this;
    }

    template<typename ELE, typename = decltype(std::declval<T &>() *= std::declval<const ELE &>())>
    SparseVec operator *(const ELE &element) const {
        SparseVec<T> result(*this);
        result *= element;
        return result;
    }

    template<typename OTHER, typename = decltype(std::declval<T &>() += std::declval<T>() * (std::declval<OTHER>()[0]))>
    T operator *(const OTHER &other) const {
        T dotProd = 0;
        for(int i=0; i < sparseSize(); ++i) {
            dotProd += values[i] * other[indices[i]];
        }
        return dotProd;
    }


    SparseVec<T> operator -() const {
        SparseVec<T> negation(sparseSize());
        for(int i=0; i<sparseSize(); ++i) {
            negation.indices[i] = indices[i];
            negation.values[i] = -values[i];
        }
        return negation;
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

private:
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & indices & values;
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


// Adds element,
// doesn't delete if ind already exists
template<class T>
void SparseVec<T>::insert(int i, const T &v) {
    if (v != 0) {
        indices.push_back(i);
        values.push_back(v);
    }
}


template<class T>
void SparseVec<T>::insert(int i, T &&v) {
        if (v != 0) {
            indices.push_back(i);
            values.push_back(std::move(v));
        }
}


template<class T>
void SparseVec<T>::insert(const SparseVec<T> &other) {
    for(int nzi=0; nzi<other.sparseSize(); ++nzi) {
        insert(other.indices[nzi], other.values[nzi]);
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
    out << "{";
    for(int i=0; i<sVector.sparseSize(); ++i) {
        out << "X[" << sVector.indices[i] << "]=" << sVector.values[i] << "  ";
    }
    out << "}";
    return out;
}


template<typename T, typename OTHER, typename = decltype(std::declval<OTHER>()[0] += std::declval<T>())>
OTHER &operator +=(OTHER &lhs, const SparseVec<T> &rhs) {
    for(int nz=0; nz < rhs.sparseSize(); ++nz) {
        lhs[rhs.indices[nz]] += rhs.values[nz];
    }
    return lhs;
}


template<typename T, typename OTHER, typename = decltype(std::declval<OTHER>()[0] -= std::declval<T>())>
OTHER &operator -=(OTHER &lhs, const SparseVec<T> &rhs) {
    for(int nz=0; nz < rhs.sparseSize(); ++nz) {
        lhs[rhs.indices[nz]] -= rhs.values[nz];
    }
    return lhs;
}


template<typename T, typename OTHER, typename = decltype(std::declval<SparseVec<T>>() * std::declval<OTHER>())>
auto operator *(const OTHER &lhs, const SparseVec<T> &rhs) {
    return rhs * lhs;
}


#endif //GLPKTEST_SPARSEVEC_H
