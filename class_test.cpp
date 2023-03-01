#include <iostream>

#include "class_template.h"

using namespace std;
  
int main()
{
    int arr[5] = { 1, 2, 3, 4, 5 };
    Array a(arr, 5);
    a.print();
    long double arr2[5] = { 0.1, 0.2, 0.3, 0.4, 0.6 };
    Array a2(arr2, 5);
    a2.print();
    return 0;
}