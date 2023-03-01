//https://www.youtube.com/watch?v=jmo_HN_-PxI
//
#ifndef QBMATRIX2_h
#define QBMATRIX2_h


// loading standard libraries
#include <cstdlib>
#include <string>

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <exception>
#include "qbVector.h" // class for vetor operations in qbLinAlg linear algebra library

template <class T>
class qbMatrix2{
    public:
    // define the various constructors
    qbMatrix2();
    qbMatrix2(int nRows,int nCols); // constructor for creating an empty matrix (all elements set to 0)
    qbMatrix2(int nRows,int nCols, const T *inputData);  // inputData is a linear array that contains the data we want to store
    qbMatrix2(int nRows,int nCols, auto inputData);  // inputData is a linear array that contains the data we want to store
    qbMatrix2(const qbMatrix2<T>& inputMatrix); // copy constructor
    qbMatrix2(int nRows, int nCols, const std::vector<T> &inputData); // constructpr for accepting vectors as input
    qbMatrix2(int nRows, int nCols, const qbVector<T> &inputData); // constructpr for accepting vectors as input

    // define the destructors
    ~qbMatrix2();// this way we free llocated memory are the code progresses

    // Configuraion methods
    bool Resize(int numRows,int numCols); // if the process is successfull it returens True
    void SetToIdentity(); //_ // Function to convert the existing matrix into an identity matrix.
    void SetToZero(); //_ // Function to convert the existing matrix into an identity matrix.

    // Element access methods
    T GetElement(int row, int col) const {return (Sub2Ind(row,col)>=0) ? m_matrixData[Sub2Ind(row,col)] : throw std::invalid_argument("bad row or column.");};
    //T GetElement(int row, int col){return ((Sub2Ind(row,col)>=0)&&(Sub2Ind(row,col)<m_nElements)) ? m_matrixData[Sub2Ind(row,col)] : 0.0;};
    //T GetElement(int row, int col){ return m_matrixData[Sub2Ind(row,col)];};
    void SetElement(int row, int col, T elementValue){(Sub2Ind(row,col)>=0) ? m_matrixData[Sub2Ind(row,col)] = elementValue : throw std::invalid_argument("bad row or column.");};    
    //void SetElement(int row, int col, T elementValue){((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0)) ? m_matrixData[Sub2Ind(row,col)] = elementValue : m_matrixData[Sub2Ind(row,col)] = 0.0;};    
    //void SetElement(int row, int col, T elementValue){((Sub2Ind(row,col)>=0)&&(Sub2Ind(row,col)<m_nElements)) ? m_matrixData[Sub2Ind(row,col)] = elementValue : m_matrixData[Sub2Ind(row,col)] = 0.0;};    
    //void SetElement(int row, int col, T elementValue){m_matrixData[Sub2Ind(row,col)] = elementValue;};    
    int GetNumRows() const {return m_nRows;};
    int GetNumCols() const {return m_nCols;};


    // Manipulation methods.
    // vector conversion
    qbVector<T> vector();
    // Compute matrix inverse.
    bool Inverse(); //_
    // Convert to row echelon form.
    qbMatrix2<T> RowEchelon(); //_
    // Return the transpose.
    qbMatrix2<T> Transpose() const; //_
    
    // Compute determinant.
    T Determinant(); //_

    // Overload == (equal) Operator
    bool operator == (const qbMatrix2<T>& rhs); //_ // testing of equality
	bool Compare (const qbMatrix2<T>& matrix1, double tolerance); //_ // It allows to specify a tolerance for comparing matrices
    //
	// Overload the assignment operator.
	qbMatrix2<T> operator= (const qbMatrix2<T> &rhs);
    //
	// Overload the * operator.    
    //qbMatrix2<T> operator * (const qbMatrix2<T>& rhs);
    //
    // overload operators (friends, this methods are defined outside the class and they need to have accsess to private members) 
    template <class U> friend qbMatrix2<U> operator + (const qbMatrix2<U>& lhs,const qbMatrix2<U>& rhs); // matrix_A + matrix_B
    template <class U> friend qbMatrix2<U> operator + (const U& lhs,const qbMatrix2<U>& rhs);// scalar  + matrix_B
    template <class U> friend qbMatrix2<U> operator + (const qbMatrix2<U>& lhs,const U& rhs);// matrix_A + scalar
    //
    template <class U> friend qbMatrix2<U> operator - (const qbMatrix2<U>& lhs,const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator - (const U& lhs,const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator - (const qbMatrix2<U>& lhs,const U& rhs);
    //
    template <class U> friend qbMatrix2<U> operator * (const qbMatrix2<U>& lhs,const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator * (const U& lhs,const qbMatrix2<U>& rhs);
    template <class U> friend qbMatrix2<U> operator * (const qbMatrix2<U>& lhs,const U& rhs);
    //
    // qbMatrix2 * qbVector.
    template <class U> friend qbVector<U> operator* (const qbMatrix2<U>& lhs, const qbVector<U>& rhs); //

	bool Separate(qbMatrix2<T> &matrix1, qbMatrix2<T> &matrix2, int colNum); //
	bool Join(const qbMatrix2<T>& matrix2); //
	qbMatrix2<T> FindSubMatrix(int rowNum, int colNum); //
	
	// Function to return the rank of the matrix.
	int Rank();		//
	
	bool IsSquare();    //
	bool IsRowEchelon();	//
	bool IsNonZero();	//
	bool IsSymmetric(); //
	void PrintMatrix(); //
	void PrintMatrix(int precision);    //
    //
    private:
        int Sub2Ind(int row, int col) const {return ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0)) ? ((row*m_nCols)+col) : -1;};
		bool CloseEnough(T f1, T f2);   //
		void SwapRow(int i, int j); //
		void MultAdd(int i, int j, T multFactor); //
		void MultRow(int i, T multFactor);
		int FindRowWithMaxElement(int colNumber, int startingRow);	
    //
    private:
        T *m_matrixData;
        int  m_nRows, m_nCols, m_nElements;
};

// Constructor that takes a linear vector as input
template <class T> 
qbMatrix2<T>::qbMatrix2(int nRows,int nCols, auto inputData)
{
// check that the number of input elements is consistent with the size of the matrix
m_nRows = nRows;
m_nCols = nCols;
m_nElements = m_nRows*m_nCols;
m_matrixData = new T[m_nElements];
for(int ii = 0;ii<m_nElements;ii++){
    m_matrixData[ii] = inputData[ii];
    //cout << m_matrixData[ii]<< " ";
    }
//
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////// Operator Overloading //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// matrix_A + matrix_B
template <class T> 
qbMatrix2<T> operator + (const qbMatrix2<T>& lhs,const qbMatrix2<T>& rhs)
{
int numRows = lhs.m_nRows;
int numCols = lhs.m_nCols;
int numElemetns = lhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;ii++){
        tempResult[ii] = lhs.m_matrixData[ii] + rhs.m_matrixData[ii]; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}

// scalar + matrix_B
template <class T> 
qbMatrix2<T> operator + (const T& lhs,const qbMatrix2<T>& rhs)
{
int numRows = rhs.m_nRows;
int numCols = rhs.m_nCols;
int numElemetns = rhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs + rhs.m_matrixData[ii]; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}

// matrix_A + scalar
template <class T> 
qbMatrix2<T> operator + (const qbMatrix2<T>& lhs,const T& rhs)
{
int numRows = lhs.m_nRows;
int numCols = lhs.m_nCols;
int numElemetns = lhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs.m_matrixData[ii] + rhs; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}

// matrix_A - matrix_B
template <class T> 
qbMatrix2<T> operator - (const qbMatrix2<T>& lhs,const qbMatrix2<T>& rhs)
{
int numRows = lhs.m_nRows;
int numCols = lhs.m_nCols;
int numElemetns = lhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;ii++){
        tempResult[ii] = lhs.m_matrixData[ii] - rhs.m_matrixData[ii]; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}

// scalar - matrix_B
template <class T> 
qbMatrix2<T> operator - (const T& lhs,const qbMatrix2<T>& rhs)
{
int numRows = rhs.m_nRows;
int numCols = rhs.m_nCols;
int numElemetns = rhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs - rhs.m_matrixData[ii]; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}




// matrix_A - scalar
template <class T> 
qbMatrix2<T> operator - (const qbMatrix2<T>& lhs,const T& rhs)
{
int numRows = lhs.m_nRows;
int numCols = lhs.m_nCols;
int numElemetns = lhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs.m_matrixData[ii] - rhs; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}


// matrix * vector
template <class T>
qbVector<T> operator* (const qbMatrix2<T>& lhs, const qbVector<T>& rhs)
{
	// Verify the dimensions of the inputs.
	if (lhs.m_nCols != rhs.GetNumDims())
		throw std::invalid_argument("Number of columns in matrix must equal number of rows in vector.");
	
	// Setup the vector for the output.
	qbVector<T> result(lhs.m_nRows);
	
	// Loop over rows and columns and perform the multiplication operation element-by-element.
	for (int row=0; row<lhs.m_nRows; ++row)
	{
		T cumulativeSum = static_cast<T>(0.0);
		for (int col=0; col<lhs.m_nCols; ++col)
		{
			cumulativeSum += (lhs.GetElement(row,col) * rhs.GetElement(col));
		}
		result.SetElement(row, cumulativeSum);
	}
	
	return result;
}


// matrix_A * matrix_B
template <class T> 
qbMatrix2<T> operator * (const qbMatrix2<T>& lhs,const qbMatrix2<T>& rhs)
{
int l_numRows = lhs.m_nRows;
int l_numCols = lhs.m_nCols;
int r_numRows = rhs.m_nRows;
int r_numCols = rhs.m_nCols;
//int numElemetns = rhs.m_nElements;
int numElemetns = l_numRows*r_numCols;
//
T elementResult = 0.0;
int lhsKK;
int rhsKK;
int nn;
//
//
if(l_numCols==r_numRows){
    // standar matrix multiplication contition
    // the results will have the same number or rows as the lhs and the same number of columns as the rhs
    T *tempResult = new T[numElemetns];
    //
    // loop though each row of the lhs
    for(int ii=0;ii<l_numRows;ii++){
        //
        // loop though each column of the rhs
        for(int jj=0;jj<r_numCols;jj++){
            //
            elementResult = 0.0;
            // loop though each element of the rhs row and lhs columns
            for(int kk=0;kk<r_numRows;kk++){
                //lhsKK = lhs.Sub2Ind(ii,kk);// linear index of lhs
                //rhsKK = rhs.Sub2Ind(kk,jj);// linear index of rhs   
                lhsKK = (ii*l_numCols)+kk;
                rhsKK = (kk*r_numCols)+jj;
                elementResult += ( lhs.m_matrixData[lhsKK]*rhs.m_matrixData[rhsKK] );
            }
            nn = ii*r_numCols+jj;//result Linear Index
            tempResult[nn] = elementResult;
        }        
    }
    qbMatrix2<T> result(l_numRows,r_numCols,tempResult);
    delete[] tempResult;
    return result;  
 }
 else{
     qbMatrix2<T> result(1,1);
     return result;  
 }

}//// matrix_A * matrix_B



// scalar * matrix_B
template <class T> 
qbMatrix2<T> operator * (const T& lhs,const qbMatrix2<T>& rhs)
{
int numRows = rhs.m_nRows;
int numCols = rhs.m_nCols;
int numElemetns = rhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs * rhs.m_matrixData[ii]; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}

// matrix_A * scalar 
template <class T> 
qbMatrix2<T> operator * (const qbMatrix2<T>& lhs,const T& rhs)
{
int numRows = lhs.m_nRows;
int numCols = lhs.m_nCols;
int numElemetns = lhs.m_nElements;
T *tempResult = new T[numElemetns];
for(int ii=0;ii<numElemetns;++ii){
        tempResult[ii] = lhs.m_matrixData[ii] * rhs; 
    }
    qbMatrix2<T> result(numRows,numCols,tempResult);
    delete[] tempResult;
    return result;  
}


#endif 

