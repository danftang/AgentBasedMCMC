//
// Created by daniel on 04/12/2021.
//

#include "MutableCategoricalArray.h"

std::uniform_real_distribution<double> MutableCategoricalArray::uniformDist(0.0,1.0);

void MutableCategoricalArray::set(int index, double probability) {
    double sum = probability;
    int indexOffset = 1;
    while((indexOffset & index) == 0 && indexOffset < size()) {
        int descendantIndex = index + indexOffset;
        if(descendantIndex < size()) sum += tree[descendantIndex];
        indexOffset = indexOffset << 1;
    }
    double delta = sum - tree[index];
    int ancestorIndex = index;
    tree[index] = sum;
    while(indexOffset < size()) {
        ancestorIndex = ancestorIndex ^ indexOffset;
        tree[ancestorIndex] += delta;
        do {
            indexOffset = indexOffset << 1;
        } while((ancestorIndex & indexOffset) == 0 && indexOffset < size());
    }
}

double MutableCategoricalArray::descendantSum(int index) const {
    int indexOffset = 1;
    double sum = 0.0;
    while((indexOffset & index) == 0 && indexOffset < size()) {
        int descendantIndex = index + indexOffset;
        if(descendantIndex < size()) sum += tree[descendantIndex];
        indexOffset = indexOffset << 1;
    }
    return sum;
}

int MutableCategoricalArray::highestOneBit(int i) {
    i = i | (i >> 1);
    i = i | (i >> 2);
    i = i | (i >> 4);
    i = i | (i >> 8);
    i = i | (i >> 16);
    return i - (i >> 1);
}

