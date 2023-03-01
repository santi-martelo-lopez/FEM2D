#ifndef CLASS_TEMPLATE_H
#define CLASS_TEMPLATE_H

#include <iostream>
using namespace std;
  
template <typename T> class Array {
private:
    T* ptr;
    int size;
    T aa = 0.1;
  
public:
    Array(T arr[], int s);
    void print();
};
#endif