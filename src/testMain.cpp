//
// Created by daniel on 18/03/2022.
//

#include <map>
#include <iostream>
#include <functional>

class MyClass {
public:
    int x;

    MyClass() { std::cout << "MyClass default constructor" << std::endl;}
    MyClass(const MyClass &copy) { std::cout << "MyClass copy constructor" << std::endl;}
    MyClass(MyClass &&move) { std::cout << "MyClass move constructor" << std::endl;}

    int myFunc() { return x; }

};

class MyDerived: public MyClass {
};

bool myFunc(const MyClass &m) {
    return m.x < 4;
}

MyDerived myFunc2() {
    return MyDerived();
}


int main(int argc, char *argv[]) {
    std::function<const MyClass &()> f;
    f = myFunc2;
    std::cout << f().x << std::endl;
}