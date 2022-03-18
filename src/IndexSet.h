// Represents a subset of [0,N].
// Allows O(1) member insertion, O(1) membership testing and efficient iteration
// but no deletion other than clear()
//
// Works by having a vector of bools and a vector of members, so N shouldn't be
// too big.
//
// Created by daniel on 17/03/2022.
//

#ifndef ABMCMC_INDEXSET_H
#define ABMCMC_INDEXSET_H

#include <vector>

class IndexSet {
protected:
    std::vector<bool>   _isMember;  // is a given integer a member?
    std::vector<int>    _members;   // list of members

public:
    explicit IndexSet(size_t maxIndex) : _isMember(maxIndex+1, false) { }

    void insert(int index) {
        if(!_isMember[index]) {
            _isMember[index] = true;
            _members.push_back(index);
        }
    }

    void clear() {
        for(int index : _members) _isMember[index] = false;
        _members.clear();
    }

    size_t size() { return _members.size(); } // number of members
    size_t maxIndex() { return _isMember.size()-1; } // number of members

    auto begin() { return _members.begin(); }
    auto end() { return _members.end(); }
    auto rbegin() { return _members.rbegin(); }
    auto rend() { return _members.rend(); }

    bool isMember(int index) { return _isMember[index]; }

};

#endif //ABMCMC_INDEXSET_H
