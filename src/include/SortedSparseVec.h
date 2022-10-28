//
// Created by daniel on 21/10/22.
//

#ifndef ABMCMC_SORTEDSPARSEVEC_H
#define ABMCMC_SORTEDSPARSEVEC_H

#include <vector>
#include <map>
#include <algorithm>

template<class T>
class SortedSparseVec : public std::vector<std::pair<int, T>> {
public:
    typedef std::vector<std::pair<int,T>> vector_type;
    using typename vector_type::iterator;
    using typename vector_type::const_iterator;
    using typename vector_type::value_type;
    using vector_type::erase;

    static const T zero;

    SortedSparseVec(): std::vector<std::pair<int,T>>() {}
    SortedSparseVec(const std::map<int, T> &map) : std::vector<std::pair<int, T>>(map.begin(), map.end()) {}

    class MutableAccessor {
    public:
        SortedSparseVec<T> &sparseVec;
        const int index;

        MutableAccessor(SortedSparseVec<T> &sparseVec, int index): sparseVec(sparseVec), index(index) {}
        inline operator const T &() const { return sparseVec.get(index); }
        inline const T &operator =(const T &value) const { sparseVec.set(index, value); return value; }
    };

    inline MutableAccessor operator [](int index) { return MutableAccessor(*this, index); }
    inline const T & operator [](int index) const { return get(index); }


    inline void set(int index, const T &value) {
        iterator lb = find(index);
        if(value != 0) {
            if (lb != this->end()) {
                if (lb->first == index) {
                    lb->second = value;
                } else {
                    lb = this->emplace(lb, index, value);
                }
            } else {
                this->emplace_back(index, value);
            }
        } else if(lb != this->end() && lb->first == index) {
            vector_type::erase(lb);
        }
    }


    inline const T &get(int index) const {
        const_iterator lb = find(index);
        return(lb != this->end() && lb->first == index) ? lb->second:zero;
    }

//    bool erase(int index) {
//        iterator lb = find(index);
//        if(lb->first != index) return false;
//        vector_type::erase(lb);
//        return true;
//    }

    iterator find(int index) {
        return std::lower_bound(this->begin(), this->end(), index, [](const std::pair<int,T> &entry, int index) {
            return entry.first < index;
        });
    }

    const_iterator find(int index) const {
        return std::lower_bound(this->begin(), this->end(), index, [](const std::pair<int,T> &entry, int index) {
            return entry.first < index;
        });
    }

    SortedSparseVec &operator *=(const T &constant) {
        if(constant != 0) {
            for (value_type &entry: *this) entry.second *= constant;
        } else {
            this->clear();
        }
        return *this;
    }

    SortedSparseVec &operator /=(const T &constant) {
        assert(constant != 0);
        for(value_type &entry: *this) entry.second /= constant;
        return *this;
    }

    SortedSparseVec &operator+=(const SortedSparseVec<T> &other) {
        if(other.empty()) return *this;
        std::vector<std::pair<int, T>> result;
        result.reserve(this->size() + other.size());
        auto extractionPointA = this->begin();
        auto extractionPointB = other.begin();
        while (extractionPointA != this->end() && extractionPointB != other.end()) {
            const auto &x = *extractionPointA;
            const auto &y = *extractionPointB;
            if (x.first < y.first) {
                result.push_back(x);
                ++extractionPointA;
            } else if (y.first < x.first) {
                result.push_back(y);
                ++extractionPointB;
            } else {
                T sum = x.second + y.second;
                if(sum != 0) result.template emplace_back(x.first, sum);
                ++extractionPointA;
                ++extractionPointB;
            }
        }
        while (extractionPointA != this->end()) result.push_back(*extractionPointA++);
        while (extractionPointB != other.end()) result.push_back(*extractionPointB++);
        std::vector<std::pair<int, T>>::operator=(std::move(result));
        return *this;
    }

    SortedSparseVec &weightedPlusAssign(const T &weight, const SortedSparseVec<T> &other) {
        if(weight == 0 || other.empty()) return *this;
        std::vector<std::pair<int, T>> result;
        result.reserve(this->size() + other.size());
        auto extractionPointA = this->begin();
        auto extractionPointB = other.begin();
        while (extractionPointA != this->end() && extractionPointB != other.end()) {
            const auto &x = *extractionPointA;
            const auto &y = *extractionPointB;
            if (x.first < y.first) {
                result.push_back(x);
                ++extractionPointA;
            } else if (y.first < x.first) {
                result.emplace_back(y.first, weight * y.second);
                ++extractionPointB;
            } else {
                T sum = x.second + weight * y.second;
                if(sum != 0) result.template emplace_back(x.first, sum);
                ++extractionPointA;
                ++extractionPointB;
            }
        }
        while (extractionPointA != this->end()) result.push_back(*extractionPointA++);
        while (extractionPointB != other.end()) {
            const std::pair<int,T> & y = *(extractionPointB++);
            result.emplace_back(y.first, weight * y.second);
        }
        std::cout << *this << " + " << weight << " * " << other << " = " << result << std::endl;
        std::vector<std::pair<int, T>>::operator=(std::move(result));
        return *this;
    }

    SortedSparseVec &operator =(const SparseVec<T> &unsortedSparseVec) {
        for(const auto &entry: unsortedSparseVec) {
            this->template emplace_back(entry.index(), entry.value());
        }
        sort();
        return *this;
    }

    void sort() {
        std::sort(this->begin(), this->end(), [](const std::pair<int,T> &a, const std::pair<int,T> &b) {
            return a.first < b.first;
        });
    }


};

template<class T> const T SortedSparseVec<T>::zero = 0;

#endif //ABMCMC_SORTEDSPARSEVEC_H
