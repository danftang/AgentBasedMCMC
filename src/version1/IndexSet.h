// Represents a subset of [0,N].
// Allows O(1) member insertion, O(1) membership testing and efficient iteration
// but no deletion other than clear()
//
// Also maintains the order of insertion, and provides an operator[] to access
// the i'th inserted member.
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
    explicit IndexSet(size_t domainSize) : _isMember(domainSize, false) { }

    // returns true if the index was inserted, false if the index was already a member
    bool insert(int index) {
        if(!_isMember[index]) {
            _isMember[index] = true;
            _members.push_back(index);
            return true;
        }
        return false;
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

    int operator[](int indexOfMember) { return _members[indexOfMember]; }

};

#endif //ABMCMC_INDEXSET_H
