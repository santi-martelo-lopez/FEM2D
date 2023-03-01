#ifndef QUAD_H
#define QUAD_H

#include <math.h>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>

#include <fstream> // reading / writing files
#include <iomanip> // for controlling output text alingment with setw(nn)

#include <exception>

#include "qbVector.h" // linear algebra library
#include "qbmatrix2.h" // linear algebra library
//#include "qbLinSolve.h" // linear algebra library

#include <cblas.h>
#include <lapacke.h>
//#include "cblas.h"

///////////////////////////// Inverse Matrix
#include <bits/stdc++.h>
#define ND 2     //Number of coordinate dimensions
#define NF 4    //Number of shape functions
#define NN 4    //Number of nodes per element
#define NEV 8    //nev = nn*nd, Number of variables in an element
#define dNF 8   // Number of shape function derivatives
#define NS 3    // Number of stress components
#define NG 4    // Number of Gaussian Integration Points
/////////////////////////////

using namespace std;

template <class T> 
class quad
{
public:
// constructor
quad();
quad(T input_eta,T input_xi);
//quad(T input_eta,T input_xi, const T *inputCoord);
quad(T input_eta,T input_xi, T inputX[NN], T inputY[NN]);
/////////////////////////////////    // Element access methods
T get_eta() {return this->eta;};
T get_xi() {return this->xi;};
//
// shape function values
vector<T> get_sF_vec() const {return sF_vec;}; // 1x4 matrix in vector form
vector<T> get_local_sf_vec() const {return local_sf_vec;}; //2x8 matrix in vector form
//
// shape function derivatives
vector<T> get_dsF_vec() const {return dsF_vec;}; // 2x4 matrix in vector form
vector<T> get_dsFdeta_vec() const {return dsFdeta_vec;}; // 1x4 matrix in vector form
vector<T> get_dsFdxi_vec() const {return dsFdxi_vec;}; // 1x4 matrix in vector form
// Jacobian matrix
vector<T> get_Jm2D_vec()const {return (J_flag == 1) ? Jm2D_vec :  throw std::invalid_argument("Singular matrix.");};
// Inverse Jacobian matrix
vector<T> get_invJm2D_vec()const {return (J_flag == 1) ? invJm2D_vec : throw std::invalid_argument("Singular matrix.");};
// Cartesian Derivatives
vector<T> get_global_sd_vec()const {return (J_flag == 1) ? global_sd_vec : throw std::invalid_argument("Singular matrix.");};
// B matrix
vector<T> get_B_vec() const {return (J_flag == 1) ? B_vec : throw std::invalid_argument("Singular matrix.");};
vector<T> get_coord_vec() const {return coord_vec;};
vector<T> get_u_vec() const {return u_vec;};
// D matrix
vector<T> get_D_vec() const {return D_vec;};
T get_pr1() const {return pr1;};
T get_pr2() const {return pr2;};
T get_alpha() const {return alpha;};
/////////////////////////////////
void Eval_local_sf(T eta, T xi);
void Eval_local_sd(T eta, T xi);
void Eval_local_sd();
void Eval_Jm2D(T x[NN], T y[NN]);
void Eval_Jm2D(vector<T> xy);
void Eval_Bm2D();
void Eval_inv_Jm2D(){(!inverse(Jm2D, invJm2D) == true) ? J_flag=0 : J_flag=1;}
void Eval_global_sd();
void Eval_D_vec(T E, T nu);
void print_sf();
void print_local_sd();
void print_Jm2D();
void print_invJm2D();
void print_global_sd();
void print_bm2D();
void Eval_principalStress(qbMatrix2<T> estress); // computes principal stresses and angle. estress = [sxx, syy, tauxy]
///////////////////////////////// Inverse Matrix
void getCofactor(T A[ND][ND], T temp[ND][ND], int p, int q,int n);
//
/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
T determinant(T A[ND][ND], T n);
//s
void adjoint(T A[ND][ND], T adj[ND][ND]);
bool inverse(T A[ND][ND], T inverse[ND][ND]);
/////////////////////////////////

//deconstructor
//virtual ~quad();
~quad();

protected:

private:
const int nodes = NN; // number of nodes
T eta, xi;
vector<T> sF_vec;
vector<T> local_sf_vec;
vector<T> dsF_vec;
vector<T> dsFdeta_vec;
vector<T> dsFdxi_vec;
vector<vector<T>> local_sd_vec;   // local derivatives of the shape functions
vector<vector<T>> global_sd = { {0, 0, 0, 0},{0, 0, 0, 0}}; // Cartesian Derivatives of the shape functions
vector<T> global_sd_vec; //cartesian global derivatives of the shape functions;
T Jm2D[ND][ND] = {{0,0},{0,0}};
vector<T> Jm2D_vec;
T invJm2D[ND][ND] = {{0,0},{0,0}};
vector<T> invJm2D_vec;
int J_flag;
T bm2D[NS][dNF];
//vector<vector<T>> B_vec;
vector<T> B_vec;
vector<T> D_vec;
vector<T> coord_vec;
vector<T> u_vec;
T pr1,pr2,alpha; // principal stress anlong new x and y axis and angle of the new coordinate system
};


template <class T> 
class quadModel
{
public:
// constructor
quadModel();
quadModel(string elm_data_FileName);
quadModel(T E_val,T nu_val, T h_val);


//Modeller methods
vector<T> get_D_vec() const {return D_vec;};
vector<T> get_sF_vals(T xi, T eta);  //returns local shape function values
qbMatrix2<T> get_dsF_matrix(T xi, T eta); //returns local shape function derivative  values
qbMatrix2<T> Bmtx_quad(T xi, T eta, vector<T> xye,  T &J1Det); //returns the B matrix for a quad linear element
qbMatrix2<T> Bmtx_quad(T xi, T eta, vector<T> xye); //returns the B matrix for a quad linear element
qbMatrix2<T> Ematx_quad(int nev, int ng, vector<T> xye, qbMatrix2<T> D, T t);//function Ke=Ematx(nev,ng,xyel,D,t)
/*
nev = nn*nd : number of variables in an element
nn : number of nodes per element
nd : number of coordinate dimensions / degrees of freedom
ns : number of stress components
ng : number of Gaussian integration points
t : element thickness
*/
//
//
void Eval_D_vec(T E, T nu); // Element Stress Matrix
//
//
virtual void basicAssembler();         //returns the B amtrix and stress
virtual void basicElementAssembler();  //returns the local Ke stiffness matrix
virtual void basicKgAssembler();       //returns the global Kg stiffness matrix
virtual void basicRhsAssembler();      //returns the global RHS vector
virtual void basicKgmRHSmSolver();     //returns the global displacements ug after solving the system [Kgm]*ug = RHSm by inverting Kgm matrix   
virtual void qbLinSolveKgmRHSmSolver();//returns the global displacements ug after solving the system [Kgm]*ug = RHSm by Gauss-Seidal substitution
virtual void generalLAPCAKSolveKgmRHSmSolver(); //returns the global displacements ug after solving the system [Kgm]*ug = RHSm by using a general LAPACKE solver LAPACKE_dgesv
virtual void bandedLAPCAKSolveKgmRHSmSolver(); //returns the global displacements ug after solving the system [Kgm]*ug = RHSm by using a LAPACKE_dgbsv LAPACK solver
void vtkData(string FileName); // writes nodal stresss, strains, deformation, initial and final location
void Stresses(); // computes nodal stresss, strains, deformation
//
int qbLinSolve(const qbMatrix2<T> &aMatrix, const qbVector<T> &bVector, qbVector<T> &resultVec);
//
// xi, eta: local coodinates for interpolation
// xy: globl nodal coordinates

//deconstructor
//virtual ~quadModel();
~quadModel();



protected:

private:
// Material Properties
T E,nu; // E: Young's Modulus, nu: Poisson Ration, h: element thickness
T mu;
vector<T> D_vec;
//
// element map
// check unsigned long int
long int ne = 0; //  total number of elements
long int nds = 0; // total number of nodes
long int nev = 0; // total number of degrees of freedom in x and y axis
long int nevf; // total number of degrees of freedom in x, y and xy axis
long int nf = 0; // number of fixed nodes
long int ndsd = 0; // number of dynamic constrains
int selectFlag;
int modelType; // type of element, 4: quad linear
T gt; // element thicknes
qbMatrix2<int> elementMap; // global element map of nodes
qbMatrix2<int> nofElementMap; // global nodal indexes map where constrains are imposed
qbMatrix2<T> ifp; // global nodal map of constrains values: 1: fixed node, 0: free node
qbMatrix2<int> fxyElementMap; // global nodal indexes map where forces are imposed
qbMatrix2<T> pre; // prescribed displacements
qbMatrix2<T> fxy,RHSm; // global nodal map of force values
//qbVector<T> RHSmVec; // global nodal map of force values
qbMatrix2<T> Kg,Kgm;  // global stiffness matrix

//
// node coordinates
qbMatrix2<T> gXYvec,gpXYvec;     //global nodal coordinates, global coordinates of gauss points
vector <T> eXYvec1,eXYvec2;     // test local coordinates 
// node displacements
qbMatrix2<T> gu,gstrain,gstress,gpstress;        // global nodal displacements, strains, stresses and princial stresses 
vector <T> ue_xy1,ue_xy2;  // test local displacements
};

#endif