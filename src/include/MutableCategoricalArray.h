// This class represents a categorical probability distribution over an integer range 0..N
// where each integer is associated with a weight, w_i, and a probability given by
//
// P(i) = w_i / \sum_j w_j
//
// The weights can be read and modified using the bracket operator [],
// or alternatively by using the get() and set() member functions. Modification runs in
// O(log(N)) time while reading takes amortized O(1) time and worst case O(log(N)) time.
// The size of the array can be increased or decreased in O(log(N)) time with push_back()
// and pop_back().
//
// A random draw from the distribution can be taken using the call operator () with
// a random number generator (e.g. std::mt19937). This also runs in O(log(N)) time.
//
// The sum of all weights can be accessed in O(1) time using sum()
//
// If all probabilities need modifying simultaneously, this can be done in O(N) time using
// the setAll method.
//
// Internally this is stored as a binary sum tree. However, we
// only store sums for the root node and nodes that are right-hand children. This allows
// us to store the tree in an array of doubles of the same size as the number of leaf nodes
// while still uniquely identifying the complete binary tree.
// This is done by noting that the binary tree has leaves that are all at the same level, and
// that the leaf nodes can be mapped to the integers 0...N by choosing the integer whose binary
// digits correspond to the sequence of left/right (0/1) branches on the path from root to leaf.
// Each node that is a right child is then associated with the leaf that is reached by taking
// only left children from that node. So, given an array, each right child
// in a sum tree can be mapped to an array entry whose index can be calculated,
// starting from the most significant bit, by following the sequence
// of left/right (0/1) branches on the path from the root to that node, then taking only
// left children to the associated leaf. Since the path to a node must end with a right branch (1)
// each array entry also maps to a unique right child, so we have a unique one-one mapping.
//
// This encoding allows arrays of any size (not just integer multiples of 2)
// and allows very efficient modification of probabilities and sampling in O(log(N)) time.
#ifndef CPP_MUTABLECATEGORICALARRAY_H
#define CPP_MUTABLECATEGORICALARRAY_H

#include <functional>
#include <array>
#include <random>
#include <ostream>

class MutableCategoricalArray {

    // This class allows array operator [] syntax for both reading and writing
    class EntryRef {
        int i;
        MutableCategoricalArray &p;
    public:

        EntryRef(int index, MutableCategoricalArray &dist): i(index), p(dist) { }
        operator double() const { return p.get(i); }
        double weight() const { return p.get(i); }
        double operator =(double weight) { p.set(i, weight); return weight; }
        double operator =(const EntryRef &otherRef) {   // reference assignment semantics
            double w_i = otherRef.weight();
            p.set(i, w_i); return w_i; }
    };

    std::vector<double> tree;
    int indexHighestBit;         // 2^(number of bits necessary to hold the highest index in tree).

public:

    MutableCategoricalArray(): indexHighestBit(0) { }

    MutableCategoricalArray(int size): tree(size,0.0) {
        indexHighestBit = highestOneBit(size-1);
    }

    MutableCategoricalArray(int size, std::function<double(int)> init): MutableCategoricalArray(size) {
        for(int i=size-1; i>=0; --i) tree[i] = descendantSum(i) + init(i);
    }

    MutableCategoricalArray(std::initializer_list<double> values): MutableCategoricalArray(values.size()) {
        setAll(values);
    }

    template<typename ITERATOR,
            typename std::enable_if<
                    std::is_convertible<
                            typename std::iterator_traits<ITERATOR>::iterator_category,
                            std::input_iterator_tag
                    >::value
                            >::type>
    MutableCategoricalArray(ITERATOR begin, ITERATOR end): tree(begin, end) {
        for(int i=tree.size()-1; i>=0; --i) tree[i] += descendantSum(i);
    }


    size_t size() const { return tree.size(); }

    void reserve(size_t n) { tree.reserve(n); }

    // add a new category with index size()
    void push_back(double weight) {
        int newIndex = tree.size();
        tree.push_back(0.0);
        indexHighestBit = highestOneBit(newIndex);
        set(newIndex, weight);
    }

    // remove the highest index category.
    void pop_back() {
        set(size()-1, 0.0);
        tree.pop_back();
        indexHighestBit = highestOneBit(size()-1);
    }

    // sets the weight associated with the supplied index
    EntryRef operator [](int index) { return EntryRef(index, *this); }

    // returns the weight of the supplied index.
    double operator [](int index) const { return get(index); }

    // gets the weight associated with an index
    double get(int index) const { return tree[index] - descendantSum(index); }

    // sets the weight associated with an index
    void set(int index, double weight) {
        double sum = weight;
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


    // draws a sample from the distribution in proportion to the weights
    template<typename RNG> int operator()(RNG &generator) const;


    // Sets the un-normalised probabilities of the first N integers
    // Runs in O(N) time since descendantSum runs in amortized constant
    // time (simce average number of steps is 2 for any size of tree).
    template<typename RANDOMACCESSCONTAINER>
    void setAll(const RANDOMACCESSCONTAINER &values) {
        for(int i = values.nCategories() - 1; i >= 0; --i) {
            tree[i] = descendantSum(i) + values[i];
        }
    }

    void setAll(std::initializer_list<double> values) {
        auto it = std::rbegin(values);
        for(int i = values.size()-1; i>=0; --i) {
            tree[i] = descendantSum(i) + *it++;
        }
    }

    // the sum of all weights (doesn't need to be 1.0)
    double sum() const { return size()==0?0.0:tree[0]; }

    // Returns the normalised probability of the index'th element
    double P(int index) const { return get(index) / sum(); }

    friend std::ostream &operator <<(std::ostream &out, const MutableCategoricalArray &distribution) {
        for(int i=0; i<distribution.tree.size(); ++i) {
            out << distribution[i] << " ";
        }
        return out;
    }

protected:

    // Calculates the sum of all right children associated with a given node
    // (under left-child deletion).
    double descendantSum(int index) const {
        int indexOffset = 1;
        double sum = 0.0;
        while((indexOffset & index) == 0 && indexOffset < size()) {
            int descendantIndex = index + indexOffset;
            if(descendantIndex < size()) sum += tree[descendantIndex];
            indexOffset = indexOffset << 1;
        }
        return sum;
    }

    static int highestOneBit(int i) {
        i = i | (i >> 1);
        i = i | (i >> 2);
        i = i | (i >> 4);
        i = i | (i >> 8);
        i = i | (i >> 16);
        return i - (i >> 1);
    }
};


template<typename RNG>
int MutableCategoricalArray::operator()(RNG &generator) const {
    int index = 0;
    double target = std::uniform_real_distribution<double>(0.0, sum())(generator);
    int rightChildOffset = indexHighestBit;
    while(rightChildOffset != 0) {
        int childIndex = index+rightChildOffset;
        if(childIndex < size()) {
            if (tree[childIndex] > target) index += rightChildOffset; else target -= tree[childIndex];
        }
        rightChildOffset = rightChildOffset >> 1;
    }
    return index;
}

#endif //CPP_MUTABLECATEGORICALARRAY_H
