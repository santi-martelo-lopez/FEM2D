#include "qbmatrix2.h"

// loading stadard libraries
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include<math.h> // floor


// loading name spaces
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// CONSTRUSTORS / DESTRUCTORS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T> 
qbMatrix2<T>::qbMatrix2()
{
m_nRows = 1;
m_nCols = 1;
m_nElements = 1;
m_matrixData = new T[m_nElements];
m_matrixData[0] = 0.0;
}

// constructor for creating an empty matrix (all elements set to 0)
template <class T> 
qbMatrix2<T>::qbMatrix2(int nRows,int nCols)
{
m_nRows = nRows;
m_nCols = nCols;
m_nElements = m_nRows*m_nCols;
m_matrixData = new T[m_nElements];
for(int ii = 0;ii<m_nElements;ii++){
    m_matrixData[ii] = 0.0;
    }
//
}



// Constructor that takes a linear array as input
template <class T> 
qbMatrix2<T>::qbMatrix2(int nRows,int nCols, const T *inputData)
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

/*
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
*/


// copy constructor
template <class T> qbMatrix2<T>::qbMatrix2(const qbMatrix2<T>& inputMatrix)
{
// check that the number of input elements is consistent with the size of the matrix
m_nRows = inputMatrix.m_nRows;
m_nCols = inputMatrix.m_nCols;
m_nElements = m_nRows*m_nCols;
m_matrixData = new T[m_nElements];
for(int ii = 0;ii<m_nElements;ii++){
    m_matrixData[ii] = inputMatrix.m_matrixData[ii];
    }
//
}


// Construct from std::vector.
template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const std::vector<T> &inputData)
{
	m_nRows = nRows;
	m_nCols = nCols;
	m_nElements = m_nRows * m_nCols;
	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; ++i)
		m_matrixData[i] = inputData.at(i);
}

// Construct from qbVector<T>.
template <class T>
qbMatrix2<T>::qbMatrix2(int nRows, int nCols, const qbVector<T> &inputData)
{
	m_nRows = nRows;
	m_nCols = nCols;
	m_nElements = m_nRows * m_nCols;
	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; ++i)
		m_matrixData[i] = inputData.GetElement(i);
}

// destructor
template <class T> 
qbMatrix2<T>::~qbMatrix2()
{
    if(m_matrixData != nullptr)
        delete[] m_matrixData;// because we use the "new" command to define variable "m_matrixData"
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// CONFIGURATION METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template <class T> 
bool qbMatrix2<T>::Resize(int numRows,int numCols)
{
m_nRows = numRows;
m_nCols = numCols;
m_nElements = m_nRows*m_nCols;
delete[] m_matrixData;
m_matrixData = new T[m_nElements];
if(m_matrixData != nullptr){
        for(int ii = 0;ii<m_nElements;ii++){
            m_matrixData[ii] = 0.0;
        }
        return true;
    }
    else{
        return false;
    }


}

// Function to convert the existing matrix into an identity matrix.
template <class T>
void qbMatrix2<T>::SetToIdentity()
{
	if (!IsSquare())
		throw std::invalid_argument("Cannot form an identity matrix that is not square.");
		
	for (int row=0; row<m_nRows; ++row)
	{
		for (int col=0; col<m_nCols; ++col)
		{
			if (col == row)
				m_matrixData[Sub2Ind(row,col)] = 1.0;
			else
				m_matrixData[Sub2Ind(row,col)] = 0.0;
		}
	}
}


// Function to create a matrix of zeros.
template <class T>
void qbMatrix2<T>::SetToZero()
{
	for (int row=0; row<m_nRows; ++row)
	{
		for (int col=0; col<m_nCols; ++col)
		{
			m_matrixData[Sub2Ind(row,col)] = 0.0;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// MANIPULATION METHODS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* **************************************************************************************************
SEPARATE THE MATRIX INTO TWO PARTS, AROUND THE COLUMN NUMBER PROVIDED
(Note that the output is returned into the two qbMatrix2<T> pointers in the input argument list)
/* *************************************************************************************************/
template <class T>
bool qbMatrix2<T>::Separate(qbMatrix2<T> &matrix1, qbMatrix2<T> &matrix2, int colNum)
{
	// Compute the sizes of the new matrices.
	int numRows = m_nRows;
	int numCols1 = colNum;
	int numCols2 = m_nCols - colNum;
	
	// Resize the two matrices to the proper dimensions.
	matrix1.Resize(numRows, numCols1);
	matrix2.Resize(numRows, numCols1);
	
	// Loop over the original matrix and store data into the appropriate elements of the two
	// output matrices.
	for (int row=0; row<m_nRows; ++row)
	{
		for (int col=0; col<m_nCols; ++col)
		{
			if (col < colNum)
			{
				matrix1.SetElement(row, col, this->GetElement(row, col));
			}
			else
			{
				matrix2.SetElement(row, col-colNum, this->GetElement(row, col));
			}
		}
	}
	return true;
}

/* **************************************************************************************************
JOIN TwO MATRICES TOGETHER
/* *************************************************************************************************/
template <class T>
bool qbMatrix2<T>::Join (const qbMatrix2<T>& matrix2)
{
	// Extract the information that we need from both matrices
	int numRows1 = m_nRows;
	int numRows2 = matrix2.m_nRows;
	int numCols1 = m_nCols;
	int numCols2 = matrix2.m_nCols;
	
	// If the matrices have different numbers of rows, then this operation makes no sense.
	if (numRows1 != numRows2)
		throw std::invalid_argument("Attempt to join matrices with different numbers of rows is invalid.");
	
	// Allocate memory for the result.
	// Note that only the number of columns increases.
	T* newMatrixData = new T[numRows1*(numCols1+numCols2)];
	
	// Copy the two matrices into the new one.
	int linearIndex, resultLinearIndex;
	for (int i=0; i<numRows1; ++i)
	{
		for (int j=0; j<(numCols1+numCols2); ++j)
		{
			resultLinearIndex = (i * (numCols1+numCols2)) + j;
			
			// If j is in the left hand matrix, we get data from there...
			if (j < numCols1)
			{
				linearIndex = (i * numCols1) + j;
				newMatrixData[resultLinearIndex] = m_matrixData[linearIndex];
			}
			// Otherwise, j must be in the right hand matrix, so we get data from there...
			else
			{
				linearIndex = (i * numCols2) + (j - numCols1);
				newMatrixData[resultLinearIndex] = matrix2.m_matrixData[linearIndex];
			}
		}
	}
	
	// Update the stored data.
	m_nCols = numCols1+numCols2;
	m_nElements = m_nRows * m_nCols;
	delete[] m_matrixData;
	m_matrixData = new T[m_nElements];
	for (int i=0; i<m_nElements; ++i)
		m_matrixData[i] = newMatrixData[i];
	
	delete[] newMatrixData;
	return true;
}


/* **************************************************************************************************
COMPUTE MATRIX DETERMINANT
/* *************************************************************************************************/
template <class T>
T qbMatrix2<T>::Determinant()
{
	// Check if the matrix is square.
	if (!IsSquare())
		throw std::invalid_argument("Cannot compute the determinant of a matrix that is not square.");
		
	// If the matrix is 2x2, we can just compute the determinant directly.
	T determinant;
	if (m_nRows == 2)
	{
		determinant = (m_matrixData[0] * m_matrixData[3]) - (m_matrixData[1] * m_matrixData[2]);
	}
	else
	{
		/* Otherwise we extract the sub-matrices and then recursively call this function
			until we get to 2x2 matrices. */
			
		// We will find the sub-matrices for row 0.
		// So, loop over each column.
		T cumulativeSum = 0.0;
		T sign = 1.0;
		for (int j=0; j<m_nCols; ++j)
		{
			// And find the sub-matrix for each element.
			qbMatrix2<T> subMatrix = this->FindSubMatrix(0, j);
			
			/* Cumulatively multiply the determinant of the sub-matrix by the
				current element of this matrix and the sign variable (note the
				recursive calling of the Determinant() method). */
			cumulativeSum += this->GetElement(0, j) * subMatrix.Determinant() * sign;
			sign = -sign;
		}
		determinant = cumulativeSum;
	}	
		
	return determinant;
}

/* **************************************************************************************************
COMPUTE MATRIX INVERSE (USING GAUSS-JORDAN ELIMINATION)
/* *************************************************************************************************/
template <class T>
bool qbMatrix2<T>::Inverse()
{
	// Check if the matrix is square (we cannot compute the inverse if it isn't).
	if (!IsSquare())
		throw std::invalid_argument("Cannot compute the inverse of a matrix that is not square.");
	
	// If we get to here, the matrix is square so we can continue.
	
	// Form an identity matrix with the same dimensions as the matrix we wish to invert.
	qbMatrix2<T> identityMatrix(m_nRows, m_nCols);
	identityMatrix.SetToIdentity();
	
	// Join the identity matrix to the existing matrix.	
	int originalNumCols = m_nCols;
	Join(identityMatrix);

	// Begin the main part of the process.
	int cRow, cCol;
	int maxCount = 100;
	int count = 0;
	bool completeFlag = false;
	while ((!completeFlag) && (count < maxCount))
	{
		for (int diagIndex=0; diagIndex<m_nRows; ++diagIndex)
		{		
			// Loop over the diagonal of the matrix and ensure all diagonal elements are equal to one.
			cRow = diagIndex;
			cCol = diagIndex;

			// Find the index of the maximum element in the current column.
			int maxIndex = FindRowWithMaxElement(cCol, cRow);
		
			// If this isn't the current row, then swap.
			if (maxIndex != cRow)
			{
				//std::cout << "Swap rows " << cRow << " and " << maxIndex << std::endl;
				SwapRow(cRow, maxIndex);
			}
			// Make sure the value at (cRow,cCol) is equal to one.
			if (m_matrixData[Sub2Ind(cRow,cCol)] != 1.0)
			{
				T multFactor = 1.0 / m_matrixData[Sub2Ind(cRow,cCol)];
				MultRow(cRow, multFactor);
				//std::cout << "Multiply row " << cRow << " by " << multFactor << std::endl;
			}
		
			// Consider the column.
			for (int rowIndex=cRow+1; rowIndex<m_nRows; ++rowIndex)
			{
				// If the element is already zero, move on.
				if (!CloseEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0))
				{
					// Obtain the element to work with from the matrix diagonal.
					// As we aim to set all the diagonal elements to one, this should
					// always be valid for a matrix that can be inverted.
					int rowOneIndex = cCol;
				
					// Get the value stored at the current element.
					T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];
					
					// Get the value stored at (rowOneIndex, cCol)
					T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];
					
					// If this is equal to zero, then just move on.
					if (!CloseEnough(rowOneValue, 0.0))
					{
						// Compute the correction factor.
						// (required to reduce the element at (rowIndex, cCol) to zero).
						T correctionFactor = -(currentElementValue / rowOneValue);
	
						MultAdd(rowIndex, rowOneIndex, correctionFactor);
						
						//std::cout << "Multiply row " << rowOneIndex << " by " << correctionFactor <<
						//	" and add to row " << rowIndex << std::endl;
					}
				}
			}
		
			// Consider the row.			
			for (int colIndex=cCol+1; colIndex<originalNumCols; ++colIndex)
			{			
				// If the element is already zero, move on.
				if (!CloseEnough(m_matrixData[Sub2Ind(cRow, colIndex)], 0.0))
				{
					// Obtain the element to work with from the matrix diagonal.
					// As we aim to set all the diagonal elements to one, this should
					// always be valid for a matrix that can be inverted.
					int rowOneIndex = colIndex;
				
					// Get the value stored at the current element.
					T currentElementValue = m_matrixData[Sub2Ind(cRow, colIndex)];
					
					// Get the value stored at (rowOneIndex, colIndex)
					T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, colIndex)];
					
					// If this is equal to zero, then just move on.
					if (!CloseEnough(rowOneValue, 0.0))
					{					
						
						// Compute the correction factor.
						// (required to reduce the element at (cRow, colIndex) to zero).
						T correctionFactor = -(currentElementValue / rowOneValue);				
					
						// To make this equal to zero, we need to add -currentElementValue multiplied by
						// the row at rowOneIndex.
						MultAdd(cRow, rowOneIndex, correctionFactor);
						
						//std::cout << "Multiply row " << rowOneIndex << " by " << correctionFactor <<
						//	" and add to row " << cRow << std::endl;
					}
				}
			}
		}
		
		// Separate the result into the left and right halves.
		qbMatrix2<T> leftHalf;
		qbMatrix2<T> rightHalf;
		this->Separate(leftHalf, rightHalf, originalNumCols);
		
		// When the process is complete, the left half should be the identity matrix.
		if (leftHalf == identityMatrix)
		{
			// Set completedFlag to true to indicate that the process has completed.
			completeFlag = true;
			
			// Rebuild the matrix with just the right half, which now contains the result.			
			m_nCols = originalNumCols;
			m_nElements = m_nRows * m_nCols;
			delete[] m_matrixData;
			m_matrixData = new T[m_nElements];
			for (int i=0; i<m_nElements; ++i)
				m_matrixData[i] = rightHalf.m_matrixData[i];
		}
		
		// Increment the counter.
		count++;
	}
	
	// Return whether the process succeeded or not.
	return completeFlag;
}

/* **************************************************************************************************
COMPUTE AND RETURN THE TRANSPOSE
/* *************************************************************************************************/
template <class T>
qbMatrix2<T> qbMatrix2<T>::Transpose() const
{
	// Form the output matrix.
	// Note that we reverse the order of rows and columns, as this will be the transpose.
	qbMatrix2<T> resultMatrix(m_nCols, m_nRows);
	
	// Now loop through the elements and copy in the appropriate order.
	for (int i=0; i<m_nRows; ++i)
	{
		for (int j=0; j<m_nCols; ++j)
		{
			resultMatrix.SetElement(j, i, this->GetElement(i, j));
		}
	}
	
	return resultMatrix;
}

/* **************************************************************************************************
CONVERT TO ROW ECHELON FORM (USING GAUSSIAN ELIMINATION)
/* *************************************************************************************************/
template <class T>
qbMatrix2<T> qbMatrix2<T>::RowEchelon()
{
	/* The current matrix must have at least as many columns as rows, but note that we don't
		actually require it to be square since we assume that the user may have combined a
		square matrix with a vector. They would do this, for example, if they were trying to 
		solve a system of linear equations. */
	if (m_nCols < m_nRows)
		throw std::invalid_argument("The matrix must have at least as many columns as rows.");
	
	/* Make a copy of the matrix data before we start. We do this because the procedure below
		will make changes to the stored matrix data (it operates 'in place') and we don't want
		this behaviour. Therefore we take a copy at the beginning and then we will replace the
		modified matrix data with this copied data at the end, thus preserving the original. */
	T *tempMatrixData;
	tempMatrixData = new T[m_nRows * m_nCols];
	for (int i=0; i<(m_nRows*m_nCols); ++i)
		tempMatrixData[i] = m_matrixData[i]; 
	
	// Begin the main part of the process.
	int cRow, cCol;
	int maxCount = 100;
	int count = 0;
	bool completeFlag = false;
	while ((!completeFlag) && (count < maxCount))
	{
		for (int diagIndex=0; diagIndex<m_nRows; ++diagIndex)
		{	
			// Loop over the diagonal of the matrix and ensure all diagonal elements are equal to one.
			cRow = diagIndex;
			cCol = diagIndex;
			
			// Find the index of the maximum element in the current column.
			int maxIndex = FindRowWithMaxElement(cCol, cRow);
			
			// Now consider the column.
			// Our aim is to set all elements BELOW the diagonal to zero.
			for (int rowIndex=cRow+1; rowIndex<m_nRows; ++rowIndex)
			{
				// If the element is already zero, move on.
				if (!CloseEnough(m_matrixData[Sub2Ind(rowIndex, cCol)], 0.0))
				{
					int rowOneIndex = cCol;
				
					// Get the value stored at the current element.
					T currentElementValue = m_matrixData[Sub2Ind(rowIndex, cCol)];
					
					// Get the value stored at (rowOneIndex, cCol)
					T rowOneValue = m_matrixData[Sub2Ind(rowOneIndex, cCol)];
					
					// If this is equal to zero, then just move on.
					if (!CloseEnough(rowOneValue, 0.0))
					{
						// Compute the correction factor.
						// (required to reduce the element at (rowIndex, cCol) to zero).
						T correctionFactor = -(currentElementValue / rowOneValue);
						MultAdd(rowIndex, rowOneIndex, correctionFactor);
					}
				}
			}							
		}
		
		/* Test whether we have achieved the desired result of converting the
			matrix into row-echelon form. */
		completeFlag = this->IsRowEchelon();
		
		// Increment the counter.
		count++;
	}
	
	// Form the output matrix.
	qbMatrix2<T> outputMatrix(m_nRows, m_nCols, m_matrixData);
	
	// Restore the original matrix data from the copy.
	for (int i=0; i<(m_nRows * m_nCols); ++i)
		m_matrixData[i] = tempMatrixData[i];
	
	// Delete the copy.
	delete[] tempMatrixData;
	
	return outputMatrix;
}


/* **************************************************************************************************
COMPUTE THE RANK OF THE PROVIDED MATRIX
/* *************************************************************************************************/
template <class T>
int qbMatrix2<T>::Rank()
{
	// Convert the matrix to row-echelon form.
	qbMatrix2<T> matrixCopy = this->RowEchelon();
	
	/* If this didn't work, then we compute the rank by finding
		the largest non-zero sub-matrix with a non-zero determinant.
		
		Note that this method is slower for large matrices and therefore
		it is better to use the RowEchelon method if possible. */
	int numNonZeroRows = 0;
	if (!matrixCopy.IsRowEchelon())
	{
		// Setup a std::vector to store the sub-matrices as we go.
		std::vector<qbMatrix2<T>> subMatrixVector;
		
		// Store the current matrix into the array first.
		subMatrixVector.push_back(*this);
		
		/* Loop through the subMatrixVector until either we have tested every
			sub-matrix or the completeFlag is set. */
		bool completeFlag = false;
		int subMatrixCount = 0;
		while ((subMatrixCount < subMatrixVector.size()) && (!completeFlag))
		{	
			// Extract the currentMatrix to work with.
			qbMatrix2<T> currentMatrix = subMatrixVector[subMatrixCount];
			subMatrixCount++;
			
			// Test if this matrix is non-zero.
			if (currentMatrix.IsNonZero())
			{
				// If the determinant is non-zero, then we have our result.
				T currentMatrixDet = currentMatrix.Determinant();
				if (!CloseEnough(currentMatrixDet, 0.0))
				{
					completeFlag = true;
					numNonZeroRows = currentMatrix.GetNumRows();
				}
				else
				{
					// Extract and store each sub-matrix (if larger than 2x2).
					if ((currentMatrix.GetNumRows() > 2) && (currentMatrix.GetNumCols() > 2))
					{
						for (int i=0; i<currentMatrix.GetNumRows(); ++i)
						{
							for (int j=0; j<currentMatrix.GetNumCols(); ++j)
							{
								// Extract this sub-matrix and store.
								subMatrixVector.push_back(currentMatrix.FindSubMatrix(i,j));
							}
						}
					}
				}
			}
		}
	}
	else
	{
		/* Converting to row echelon form did work, so we can simply
			count the number of non-zero rows to get the rank. */
	
		/* If we get to here, then we can assume that the matrix is now
			in row-echelon form and we can compute the rank quite easily
			as simply the number of non-zero rows. */
			
		int nRows = matrixCopy.GetNumRows();
		int nCols = matrixCopy.GetNumCols();
			
		// Loop over each row and test whether it has at least one non-zero element.
		for (int i=0; i<nRows; ++i)
		{
			// Loop over the columns of this row.
			int colSum = 0;
			for (int j=0; j<nCols; ++j)
			{
				if (!CloseEnough(matrixCopy.GetElement(i,j), 0.0))
					colSum++;
			}
			// If colSum is greater than zero, then increment numNonZeroRows.
			if (colSum > 0)
				numNonZeroRows++;
		}
	
	}
	// The rank of the matrix is simply the number of non-zero rows.
	return numNonZeroRows;
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// OVERLOAD OPERATORS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// this function became deprecated because it does not have the closeEnough method
//template <class T> 
//bool qbMatrix2<T>::operator == (const qbMatrix2<T>& rhs)
//{
//    //checking if matrices have the same size
//    if((this->m_nRows != rhs.m_nRows)&&(this->m_nCols != rhs.m_nCols))
//        return false;
//    
//    //checking if elements are equal
//    bool flag = true;
//    for(int ii=0;ii<this->m_nElements;++ii){
//        if(this->m_matrixData[ii]!= rhs.m_matrixData[ii])
//            flag = false;
//            break;
//    }
//    return flag;
//}

// testing of equality
/* **************************************************************************************************
THE == OPERATOR
/* *************************************************************************************************/
template <class T>
bool qbMatrix2<T>::operator== (const qbMatrix2<T>& rhs)
{
	// Check if the matricies are the same size, if not return false.
	if ((this->m_nRows != rhs.m_nRows) && (this->m_nCols != rhs.m_nCols))
		return false;
		
	// Check if the elements are equal.
	bool flag = true;
	for (int i=0; i<this->m_nElements; ++i)
	{
		//if (this->m_matrixData[i] != rhs.m_matrixData[i])
		if (!CloseEnough(this->m_matrixData[i], rhs.m_matrixData[i]))
			flag = false;
	}
	return flag;
}



template <class T>
bool qbMatrix2<T>::Compare(const qbMatrix2<T>& matrix1, double tolerance)
{
	// First, check that the matrices have the same dimensions.
	int numRows1 = matrix1.m_nRows;
	int numCols1 = matrix1.m_nCols;
	if ((numRows1 != m_nRows) || (numCols1 != m_nCols))
		return false;
		
	// Loop over all the elements and compute the sum-of-squared-differences.
	double cumulativeSum = 0.0;
	for (int i=0; i<m_nElements; ++i)
	{
		T element1 = matrix1.m_matrixData[i];
		T element2 = m_matrixData[i];
		cumulativeSum += ((element1 - element2) * (element1 - element2));
	}
	double finalValue = sqrt(cumulativeSum / ((numRows1 * numCols1)-1));
	if (finalValue < tolerance)
		return true;
	else
		return false;
}


/* **************************************************************************************************
THE ASSIGNMENT (=) OPERATOR
/* *************************************************************************************************/
template <class T>
qbMatrix2<T> qbMatrix2<T>::operator= (const qbMatrix2<T> &rhs)
{
	// Make sure we're not assigning to ourself.
	if (this != &rhs)
	{
		m_nRows = rhs.m_nRows;
		m_nCols = rhs.m_nCols;
		m_nElements = rhs.m_nElements;
		
		if (m_matrixData)
			delete[] m_matrixData;
		
		m_matrixData = new T[m_nElements];
		for (int i=0; i<m_nElements; i++)
			m_matrixData[i] = rhs.m_matrixData[i];	
			
		return *this;
	}
}
/* **************************************************************************************************
Multiplication (*) OPERATOR
/* *************************************************************************************************/
/*
template <class T>
qbMatrix2<T> qbMatrix2<T>:: operator * (const qbMatrix2<T>& rhs){
int l_numRows = this->m_nRows;
int l_numCols = this->m_nCols;
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
            // loop though each element of the rhs row and lhs columns
            elementResult = 0;
            for(int kk=0;kk<r_numRows;kk++){
                //lhsKK = lhs.Sub2Ind(ii,kk);// linear index of lhs
                //rhsKK = rhs.Sub2Ind(kk,jj);// linear index of rhs   
                lhsKK = (ii*l_numCols)+kk;
                rhsKK = (kk*r_numCols)+jj;
                elementResult += ( this->m_matrixData[lhsKK]*rhs.m_matrixData[rhsKK] );
            }
            nn = ii*r_numCols+jj;//result Linear Index
            tempResult[nn] = elementResult;
        }        
    }
    //std::cout <<"Operation Ok"<<std::endl;
    qbMatrix2<T> result(l_numRows,r_numCols,tempResult);
    delete[] tempResult;
    return result;  
 }
 else{
    //std::cout <<"Bad Operation"<<std::endl;
     qbMatrix2<T> result(1,1);
     return result;  
 }
}
*/

/*
template <class T>
int qbMatrix2<T>::Sub2Ind(int row, int col) const
{
	if ((row < m_nRows) && (row >= 0) && (col < m_nCols) && (col >= 0))
		return (row * m_nCols) + col;
	else
		return -1;
}
*/

// Function to test whether the matrix is square.
template <class T>
bool qbMatrix2<T>::IsSquare()
{
	if (m_nCols == m_nRows)
		return true;
	else
		return false;
}

// Function to test whether the matrix is non-zero.
template <class T>
bool qbMatrix2<T>::IsNonZero()
{
	// Loop over every element.
	int numNonZero = 0;
	for (int i=0; i<m_nElements; ++i)
	{
		// If this element is close enough to zero, then
		// increment our numNonZero counter.
		if (!CloseEnough(m_matrixData[i], 0.0))
			numNonZero++;
	}
	
	/* If the numNonZero counter is still equal to zero, then
		the matrix must be all zeros, hence we return false. 
		Otherwise we return true. */
	return (numNonZero != 0);
}

// Function to test whether the matrix is in row-echelon form.
template <class T>
bool qbMatrix2<T>::IsRowEchelon()
{
	/* We do this by testing that the sum of all the elements in the
		lower triangular matrix is zero. */
	// Loop over each row, except the first one (which doesn't need to have any zero elements).
	T cumulativeSum = static_cast<T>(0.0);
	for (int i=1; i<m_nRows; ++i)
	{
		/* Loop over the columns that correspond to the lower triangular
			matrix according to the current row. */
		for (int j=0; j<i; ++j)
		{
			// Add this element to the cumulative sum.
			cumulativeSum += m_matrixData[Sub2Ind(i, j)];
		}
	}
	
	/* If the matrix is in row-echelon form, then cumulative sum should
		still equal zero, otherwise the matrix cannot be in row-echelon form. */
	return CloseEnough(cumulativeSum, 0.0);
}

// Function to test whether the matrix is symmetric.
template <class T>
bool qbMatrix2<T>::IsSymmetric()
{
	/* First test that the matrix is square, if it is
		not, then it cannot by symmetric. */
	if (!this->IsSquare())
		return false;
		
	// Now test for symmetry about the diagonal.
	T currentRowElement = static_cast<T>(0.0);
	T currentColElement = static_cast<T>(0.0);
	bool returnFlag = true;
	int diagIndex = 0;
	while ((diagIndex < m_nCols) && returnFlag)
	{
		int rowIndex = diagIndex + 1;
		while ((rowIndex < m_nRows) && returnFlag)
		{
			currentRowElement = this->GetElement(rowIndex, diagIndex);
			currentColElement = this->GetElement(diagIndex, rowIndex);
			
			// Compare the row and column elements.
			if (!CloseEnough(currentRowElement, currentColElement))
				returnFlag = false;
			
			// Increment row index.
			rowIndex++;
		}
		
		// Increment diagIndex.
		diagIndex++;
		
	}
	
	// Return the result.
	return returnFlag;
	
}

// Function to swap rows i and j (in place).
template <class T>
void qbMatrix2<T>::SwapRow(int i, int j)
{
	// Store a tempory copy of row i.
	T *tempRow = new T[m_nCols];
	for (int k=0; k<m_nCols; ++k)
		tempRow[k] = m_matrixData[Sub2Ind(i,k)];
		
	// Replace row i with row j.
	for (int k=0; k<m_nCols; ++k)
		m_matrixData[Sub2Ind(i,k)] = m_matrixData[Sub2Ind(j,k)];
		
	// Replace row k with the tempory copy of the original row i.
	for (int k=0; k<m_nCols; ++k)
		m_matrixData[Sub2Ind(j,k)] = tempRow[k];
	
	// Tidy up after ourselves.
	delete[] tempRow;
}

// Function to add a multiple of row j to row i (in place).
template <class T>
void qbMatrix2<T>::MultAdd(int i, int j, T multFactor)
{
	for (int k=0; k<m_nCols; ++k)
		m_matrixData[Sub2Ind(i,k)] += (m_matrixData[Sub2Ind(j,k)] * multFactor);
}

// Function to the find the row with the maximum element at the column given.
// Returns the row index.
template <class T>
int qbMatrix2<T>::FindRowWithMaxElement(int colNumber, int startingRow)
{
	T tempValue = m_matrixData[Sub2Ind(startingRow, colNumber)];
	int rowIndex = startingRow;
	for (int k=startingRow+1; k<m_nRows; ++k)
	{
		if (fabs(m_matrixData[Sub2Ind(k, colNumber)]) > fabs(tempValue))
		{
			rowIndex = k;
			tempValue = m_matrixData[Sub2Ind(k, colNumber)];
		}
	}
	return rowIndex;
}

// Function to multiply a row by the given value.
template <class T>
void qbMatrix2<T>::MultRow(int i, T multFactor)
{
	for (int k=0; k<m_nCols; ++k)
		m_matrixData[Sub2Ind(i,k)] *= multFactor;
}



// A simple function to print a matrix to stdout.
template <class T>
void qbMatrix2<T>::PrintMatrix()
{
	int nRows = this->GetNumRows();
	int nCols = this->GetNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    std::cout << std::fixed << std::setprecision(3) << this->GetElement(row, col) << "  ";
    }
	std::cout << std::endl;
	}    
}

// A simple function to print a matrix to stdout, with specified precision.
template <class T>
void qbMatrix2<T>::PrintMatrix(int precision)
{
	int nRows = this->GetNumRows();
	int nCols = this->GetNumCols();
	for (int row = 0; row<nRows; ++row)
  {
	  for (int col = 0; col<nCols; ++col)
    {
	    std::cout << std::fixed << std::setprecision(precision) << this->GetElement(row, col) << "  ";
    }
	std::cout << std::endl;
	}    
}


template <class T>
bool qbMatrix2<T>::CloseEnough(T f1, T f2)
{
    return fabs(f1-f2) < 1e-9;
}

// Function to find the sub-matrix for the given element.
template <class T>
qbMatrix2<T> qbMatrix2<T>::FindSubMatrix(int rowNum, int colNum)
{
	// Create a new matrix to store the sub-matrix.
	// Note that this is one row and one column smaller than the original.
	qbMatrix2<T> subMatrix(m_nRows-1, m_nCols-1);
	
	// Loop over the elements of the existing matrix and copy to sub-matrix as appropriate.
	int count = 0;
	for (int i=0; i<m_nRows; ++i)
	{
		for (int j=0; j<m_nCols; ++j)
		{
			// If i or j correspond to rowNum or colNum, then ignore this element.
			if ((i != rowNum) && (j != colNum))
			{
				subMatrix.m_matrixData[count] = this->GetElement(i,j);
				count++;
			}
		}
	}
	
	return subMatrix;
}

template <class T>
qbVector<T> qbMatrix2<T>::vector()
{
	qbVector<T> vec(this->m_nElements);
	for(int ii = 0; ii < this-> m_nRows; ii++){
	for(int jj = 0; jj < this-> m_nCols; jj++){
		vec.SetElement(ii*jj+jj,this->GetElement(ii,jj));
	}		
	}
	return vec;
}

template class qbMatrix2<long double>;
template class qbMatrix2<double>;
template class qbMatrix2<int>;