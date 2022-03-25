//
// Created by daniel on 18/03/2022.
//

#include <map>
#include <iostream>

class MyClass {
public:
    int x;

};

MyClass myFunc() {
    MyClass m;
    m.x = 1234;
    return m;
}

int main(int argc, char *argv[]) {
    std::cout << myFunc().x << std::endl;
}