//How to compile when using Fortran functions :
// >> gfortran -c functionsLibrary.f90	// we need to compile the Fortran files separe 
// once we execute the command above, we should have a new file functionsLibrary.o 
// >> g++ -std=c++17 -c -pthread programCpp.cpp 	//we are using the 17th bersion of the compiler and using threads
// once we execute the command above, we should have a new file programCpp.o 
// now we just have to combine them together into an executable
// >> gfortran -o test programCpp.o functionsLibrary.o
// How to run the program:
// >> ./test

//Standard way to compile :
// >> g++ -std=c++17 -pthread quad.cpp domain.cpp main.cpp -o test.exe	//we are using the 17th bersion of the compiler and using threads
// >> g++ quad.cpp domain.cpp main.cpp -o test.exe 
// >> g++ -std=c++17 qbmatrix2.cpp quad.cpp domain.cpp main.cpp -o test.exe	//we are using the 17th bersion of 
// >> g++ -std=c++17 qbVector.h qbmatrix2.cpp quad.cpp domain.cpp main.cpp -o test.exe	//we are using the 17th bersion of 
// >> g++ -std=c++17 qbVector.h qbmatrix2.cpp quad.cpp domain.cpp main.cpp -o test.exe -lm -lblas -llapack -llapacke
// >> g++ -std=c++17 qbVector.h qbmatrix2.cpp quad.cpp domain.cpp main.cpp -g -o test.exe -lm -lblas -llapack -llapacke
    // -g : option to set debugging flags

//How to run the script:
// >> ./test.exe


/*
%This code output is a numerical aproximation for a plane stresss o plane
%strain problem
*/
#include<string> // for usint string variables and related operations
#include<stdio.h>
#include<iostream>
#include<math.h>

#include <fstream> // reading / writing files
#include <iomanip> // for controlling output text alingment

#include <cctype> //
#include <random> // random numbers
#include <cstdlib>
#include <ctime>

// loading FEM libraries
#include "domain.h"     // In this class are placed attributes and methods to initialize the numerical domain
#include "quad.h"       // In this class are placed attributes and methods to use quad Finite Elements
#include "qbmatrix2.h"

// loading name sampces
using namespace std;

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
// Generic function to display the matrix.  We use it to
// display both adjoin and inverse. adjoin is integer matrix
// and inverse is a float.
/*
template <class T> void display(T A[ND][ND])
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
};*/

// function overloading, to display a matrix of object type qbMatrix2
template <class T>
void display(qbMatrix2<T> A){
    int nRows = A.GetNumRows();
    int nCols = A.GetNumCols();
    //cout <<"\n\n\n Displaying matrix:\n"<< endl;
    cout << endl;
    for(int ii=0;ii<nRows;++ii){
    for(int jj=0;jj<nCols;++jj){
        cout << A.GetElement(ii,jj)<<" ";
    }        
    cout << endl;
    }
    cout << endl;
    cout << endl;
}

// funtion for randomizing matrix elements
template <class T>
void randomMatrix(qbMatrix2<T>& A){
    int nRows = A.GetNumRows();
    int nCols = A.GetNumCols();
    T randomNumber;
    unsigned int low_dist  = 1;//min(A.m_matrixData);
    unsigned int high_dist = 100;//max(A.m_matrixData);
    //
    //
    std::srand( ( unsigned int )std::time( nullptr ) );
    for(int ii=0;ii<nRows;++ii){
    for(int jj=0;jj<nCols;++jj){
        randomNumber = (low_dist + rand() % ( high_dist - low_dist ));
        A.SetElement(ii,jj,randomNumber);
    }        
    cout << endl;
    }
}
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
int main()
{

// this file contain input data regarding the mechanical properties of the solid
string elm_data_FileName1,elm_data_FileName2,elm_data_FileName3,elm_data_FileName4;
////////////////////////
cout<< "\n\nThis code output is a numerical aproximation for a plane stresss \nor a plane strain problem.\n";

/*
//quadModel<long double> model1();
elm_data_FileName1 = "elm_data1.dat";
cout<< "\nLoading input data from file: " << elm_data_FileName1 << endl;
quadModel<long double> model1(elm_data_FileName1);
//quadModel<long double> model1(10.2,0.3,10);
// Testing basic assembler
model1.basicAssembler(); // output in log1.txt
//
//
elm_data_FileName2 = "elm_data2.dat";
cout<< "\nLoading input data from file: " << elm_data_FileName2 << endl;
quadModel<long double> model2(elm_data_FileName2);
model2.basicElementAssembler(); // output in log11.txt
//
//
elm_data_FileName3 = "elm_data3.dat";
cout<< "\nLoading input data from file: " << elm_data_FileName3 << endl;
quadModel<long double> model3(elm_data_FileName3);
model3.basicKgAssembler(); // output in log101.txt
model3.basicRhsAssembler(); // output in log101.txt
model3.basicKgmRHSmSolver(); // solves the system [Kgm]*ug = RHSm
*/
//elm_data_FileName4 = "elm_data3.dat";
elm_data_FileName4 = "elm_dataSeawall1.dat";
cout<< "\nLoading input data from file: " << elm_data_FileName4 << endl;
quadModel<long double> model4(elm_data_FileName4);
model4.basicKgAssembler(); // output in log101.txt
model4.basicRhsAssembler(); // output in log101.txt
model4.basicKgmRHSmSolver(); // solves the system [Kgm]*ug = RHSm
model4.Stresses(); // Computing streains and stresses
model4.vtkData("modelSeaWall1_vtkData.vtk"); // saving results in vtk format
/*
cout<< "\nLoading input data from file: " << "elm_data102.dat" << endl;
quadModel<double> model5("elm_data102.dat");
model5.basicKgAssembler(); // output in log101.txt
model5.basicRhsAssembler(); // output in log101.txt
model5.qbLinSolveKgmRHSmSolver(); // solves the system [Kgm]*ug = RHSm
model5.Stresses(); // Computing streains and stresses
model5.vtkData("elm_data102_5.vtk"); // saving results in vtk format

cout<< "\nLoading input data from file: " << "elm_data102.dat" << endl;
quadModel<double> model6("elm_data102.dat");
model6.basicKgAssembler(); // output in log1001.txt
model6.basicRhsAssembler(); // output in log1001.txt
model6.generalLAPCAKSolveKgmRHSmSolver(); // solves the system [Kgm]*ug = RHSm// stored in log1001.txt
model6.Stresses(); // Computing streains and stresses
model6.vtkData("elm_data102_6.vtk"); // saving results in vtk format

cout<< "\nLoading input data from file: " << "elm_data102.dat" << endl;
quadModel<double> model7("elm_data102.dat");
model7.basicKgAssembler(); // output in log101.txt
model7.basicRhsAssembler(); // output in log101.txt
model7.bandedLAPCAKSolveKgmRHSmSolver(); // solves the system [Kgm]*ug = RHSm
model7.Stresses(); // Computing streains and stresses
model7.vtkData("elm_data102_7.vtk"); // saving results in vtk format
*/
//===============================================================
//===============================================================
//         ********** END OF MAIN FUNCTION ***********
//===============================================================
//===============================================================

std::cout << "\n\n\n\n\n\nPress the Enter key to continue.\n\n";
std::cin.get();

return 0;

}