#include <iostream>
#include "class_template.h"

using namespace std;
  
  
template <typename T> Array<T>::Array(T arr[], int s)
{
    ptr = new T[s];
    size = s;
    for (int i = 0; i < size; i++)
        ptr[i] = arr[i]/10;
}
  
template <typename T> void Array<T>::print()
{
    for (int i = 0; i < size; i++)
        cout << " " << *(ptr + i);
        cout << "\n"<< aa << endl;
    cout << endl;
}

template class Array<long double>;
template class Array<int>;