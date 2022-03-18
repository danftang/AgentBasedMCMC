// A std::vector of values that keeps track of which elements have
// been changed since the last clearPerturbations()
//
// Maintains the O(1) modification of elements
// and allows O(1) query whether an element has been perturbed
// and O(n) iteration over all perturbed elements, where n is the number
// of perturbed elements.
//
// Works by having a vector of bools to track which elements have been changed
// and a vector of indexes that have been changed to allow efficient iteration over
// perturbed elements
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_PERTURBEDVECTOR_H
#define ABMCMC_PERTURBEDVECTOR_H

#include <vector>

template<class T>
class PerturbedVector {
protected:
    std::vector<T>      X;                  // the vector that this object represents
    std::vector<bool>   _isPerturbed;        // is element perturbed, by index
    std::vector<int>    _perturbedIndices;   // list of indices that have been perturbed

public:
    class ElementRef {
        PerturbedVector<T> &vec;
        int index;

        ElementRef(PerturbedVector<T> &vec, int index): vec(vec), index(index) { }

        operator T &() { return vec.X[index]; }

        inline const T &operator =(const T &value) {
            vec.perturb(index, value);
            return value;
        }
    };

    inline ElementRef operator [](int index) { return ElementRef(*this, index); }

    const T &operator [](int index) const { return X[index]; }

    void perturb(int index, const T &value) {
        X[index] = value;
        if(!_isPerturbed[index]) {
            _isPerturbed[index] = true;
            _perturbedIndices.push_back(index);
        }
    }

    void clearPerturbations() {
        for(int index : _perturbedIndices) _isPerturbed[index] = false;
        _perturbedIndices.clear();
    }

    auto begin() { return X.begin(); }
    auto end() { return X.end(); }

    const std::vector<int> &perturbedIndices() { return _perturbedIndices; }

    bool isPerturbed(int index) { return _isPerturbed[index]; }
};


#endif //ABMCMC_PERTURBEDVECTOR_H
