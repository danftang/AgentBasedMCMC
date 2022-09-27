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

    template<class INDEXIT, class VALUEIT>
    class Iterator {
    protected:
        INDEXIT  indexIt;
        VALUEIT  valueIt;
    public:
        typedef const Iterator  value_type;
        typedef value_type &    reference;
        typedef value_type *    pointer;

    public: // value interface
        typename INDEXIT::reference index() const { return(*indexIt); }
        typename VALUEIT::reference value() const { return(*valueIt); }
        bool operator !=(const Iterator &other) const { return indexIt != other.indexIt; }

    public: // iterator interface
        Iterator(const INDEXIT &indexIt, const VALUEIT &valueIt): indexIt(indexIt), valueIt(valueIt) {}

        Iterator &operator ++() { ++indexIt; ++valueIt; return *this; }
        Iterator &operator --() { --indexIt; --valueIt; return *this; }
        Iterator operator ++(int) { Iterator initVal = *this; operator ++(); return initVal; }
        Iterator operator --(int) { Iterator initVal = *this; operator --(); return initVal; }
        Iterator operator +(size_t increment) const { return Iterator(indexIt + increment, valueIt + increment); }
        Iterator operator -(size_t decrement) const { return Iterator(indexIt + decrement, valueIt + decrement); }
        reference operator *() const { return *this; }
        pointer   operator ->() const { return this; }
    };

    typedef Iterator<std::vector<int>::iterator, typename std::vector<T>::iterator> iterator;
    typedef Iterator<std::vector<int>::const_iterator, typename std::vector<T>::const_iterator> const_iterator;
    typedef Iterator<std::vector<int>::reverse_iterator, typename std::vector<T>::reverse_iterator> reverse_iterator;
    typedef Iterator<std::vector<int>::const_reverse_iterator, typename std::vector<T>::const_reverse_iterator> const_reverse_iterator;

    iterator begin() { return {indices.begin(), values.begin()}; }
    iterator end() { return {indices.end(), values.end()}; }
    const_iterator cbegin() const { return {indices.cbegin(), values.cbegin()}; }
    const_iterator cend() const { return {indices.cend(), values.cend()}; }
    const_iterator begin() const { return cbegin(); }
    const_iterator end() const { return cend(); }
    reverse_iterator rbegin() { return {indices.rbegin(), values.rbegin()}; }
    reverse_iterator rend() { return {indices.rend(), values.rend()}; }
    const_reverse_iterator crbegin() const { return {indices.crbegin(), values.crbegin()}; }
    const_reverse_iterator crend() const { return {indices.crend(), values.crend()}; }
    const_reverse_iterator rbegin() const { return crbegin(); }
    const_reverse_iterator rend() const { return crend(); }


    template<typename ELE, typename = decltype(std::declval<T &>() *= std::declval<const ELE &>())>
    SparseVec &operator *=(const ELE &element) {
        for(T &value : values) value *= element;
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
        for(const auto &entry: *this) dotProd += entry.value() * other[entry.index()];
        return dotProd;
    }


    SparseVec<T> operator -() const {
        SparseVec<T> negation;
        negation.reserve(sparseSize());
        for(const auto &entry: (*this)) {
            negation.indices.push_back(entry.index());
            negation.values.push_back(-(entry.value()));
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
    for(int index: indices) {
        if(index > max) max = index;
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
    for(const auto &entry: other) insert(entry.index(), entry.value());
}

template<class T>
void SparseVec<T>::clear() {
    indices.clear();
    values.clear();
}

// to zero-based dense array
template<class T>
std::vector<T> SparseVec<T>::toDense() const { return toDense(maxNonZeroIndex() + 1); }

template<class T>
std::vector<T> SparseVec<T>::toDense(int dimension) const {
    std::vector<T> denseVec(dimension,0);
    for (const auto &entry: *this) denseVec[entry.index()] = entry.value();
    return denseVec;
}


template<class T>
std::ostream &operator<<(std::ostream &out, const SparseVec<T> &sVector) {
    out << "{";
    for(const auto &entry: sVector) {
        out << "X[" << entry.index() << "]=" << entry.value() << "  ";
    }
    out << "}";
    return out;
}


template<typename T, typename OTHER, typename = decltype(std::declval<OTHER>()[0] += std::declval<T>())>
OTHER &operator +=(OTHER &lhs, const SparseVec<T> &rhs) {
    for(const auto &entry: rhs) lhs[entry.index()] += entry.value();
    return lhs;
}


template<typename T, typename OTHER, typename = decltype(std::declval<OTHER>()[0] -= std::declval<T>())>
OTHER &operator -=(OTHER &lhs, const SparseVec<T> &rhs) {
    for(const auto entry: rhs) lhs[entry.index()] -= entry.value();
    return lhs;
}


template<typename T, typename OTHER, typename = decltype(std::declval<SparseVec<T>>() * std::declval<OTHER>())>
auto operator *(const OTHER &lhs, const SparseVec<T> &rhs) {
    return rhs * lhs;
}


#endif //GLPKTEST_SPARSEVEC_H
