#include "quad.h"

// loading stadard libraries
#include <iostream>
#include <vector>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include<math.h> // floor

//#include "qbLinSolve.h" // linear algebra library

// loading name spaces
using namespace std;


//constructor
template <class T> quad<T>::quad()
{
    
}

template <class T> quad<T>::quad(T input_eta,T input_xi)
{
eta = input_eta;
xi = input_xi;
}
//template <class T> quad<T>::quad(T input_eta,T input_xi, const T *inputCoord)
//{
//
//eta = input_eta;
//xi = input_xi;
//}
//template <class T> quad<T>::quad(T input_eta,T input_xi, const T *inputCoord, const T *input_uvec)
//{
//
//eta = input_eta;
//xi = input_xi;    
//}
template <class T> quad<T>::quad(T input_eta, T input_xi, T inputX[NN], T inputY[NN])
{
eta = input_eta;
xi = input_xi;
Eval_local_sf(eta,xi);
Eval_local_sd(eta,xi);
Eval_Jm2D(inputX, inputY);
Eval_inv_Jm2D();
for(int ii=0;ii<ND;ii++){
for(int jj=0;jj<ND;jj++){
    invJm2D_vec.push_back(invJm2D[ii][jj]);
}        
}
Eval_global_sd();
Eval_Bm2D();
for(int ii = 0;ii<NF;ii++){
    coord_vec.push_back(inputX[ii]);
    coord_vec.push_back(inputY[ii]);
}
}

// Function for calculating the local shape functions
template <class T> void quad<T>::Eval_local_sf(T eta, T xi)
{
    T N1,N2,N3,N4;
    N1 =0.25*(1-eta)*(1-xi);
    N2 =0.25*(1+eta)*(1-xi);
    N3 =0.25*(1+eta)*(1+xi);
    N4 =0.25*(1-eta)*(1+xi);
    // Initialization
    if(sF_vec.empty()!=true){
    while (sF_vec.size() > 0) {
        sF_vec.pop_back();
    }
    }
    if(local_sf_vec.empty()!=true){
    while (local_sf_vec.size() > 0) {
            local_sf_vec.pop_back();
    }
    }
    //
    sF_vec.push_back(N1);
    sF_vec.push_back(N2);
    sF_vec.push_back(N3);
    sF_vec.push_back(N4);
    //
    for(int ii=0;ii<NF;ii++){
        local_sf_vec.push_back(sF_vec[ii]);
        local_sf_vec.push_back(0.0);
    }
    for(int ii=0;ii<NF;ii++){
        local_sf_vec.push_back(0.0);
        local_sf_vec.push_back(sF_vec[ii]);
    }
    //cout<<"   ";
    //cout<<" N1: "<<N1;
    //cout<<" N2: "<<N2;
    //cout<<" N3: "<<N3;
    //cout<<" N4: "<<N4;
}

// Function for calculating the local derivatives of the shape functions
template <class T> void quad<T>::Eval_local_sd(T eta, T xi)
{
    T dN1deta,dN2deta,dN3deta,dN4deta;
    T dN1dxi,dN2dxi,dN3dxi,dN4dxi;
    dN1deta =0.25*(-1)*(1-xi);
    dN2deta =0.25*(+1)*(1-xi);
    dN3deta =0.25*(+1)*(1+xi);
    dN4deta =0.25*(-1)*(1+xi);
    dN1dxi =0.25*(1-eta)*(-1);
    dN2dxi =0.25*(1+eta)*(-1);
    dN3dxi =0.25*(1+eta)*(+1);
    dN4dxi =0.25*(1-eta)*(+1);
    //
    // Initialization
    if(dsF_vec.empty()!=true){
    while (dsF_vec.size() > 0) {
        dsF_vec.pop_back();
    }
    }
    if(dsFdeta_vec.empty()!=true){
    while (dsFdeta_vec.size() > 0) {
        dsFdeta_vec.pop_back();
        dsFdxi_vec.pop_back();
    }
    }
    if(local_sd_vec.empty()!=true){
        for(int ii=0;ii<ND;ii++)
            local_sd_vec.pop_back();
    }
    //
    dsF_vec.push_back(dN1deta);
    dsF_vec.push_back(dN2deta);
    dsF_vec.push_back(dN3deta);
    dsF_vec.push_back(dN4deta);
    dsF_vec.push_back(dN1dxi);
    dsF_vec.push_back(dN2dxi);
    dsF_vec.push_back(dN3dxi);
    dsF_vec.push_back(dN4dxi);
    /////////////////////
    dsFdeta_vec.push_back(dN1deta);
    dsFdeta_vec.push_back(dN2deta);
    dsFdeta_vec.push_back(dN3deta);
    dsFdeta_vec.push_back(dN4deta);
    //////////////////////
    dsFdxi_vec.push_back(dN1dxi);
    dsFdxi_vec.push_back(dN2dxi);
    dsFdxi_vec.push_back(dN3dxi);
    dsFdxi_vec.push_back(dN4dxi);
    //////////////////////
    local_sd_vec.push_back(dsFdeta_vec);
    local_sd_vec.push_back(dsFdxi_vec);
    //////////////////////
}

template <class T> void quad<T>::print_sf()
//void quad::print_Eval_local_sf()
{
    cout<<" Shape Funtion values: ";
    for (auto ii=sF_vec.begin();ii!=sF_vec.end();++ii){
        cout<<" "<< *ii;
    }
   cout<<" N2 matrix: "<<endl;
    for (auto ii=0;ii<ND;++ii){
    for (auto jj=0;jj<dNF;++jj){
        cout<<" "<< local_sf_vec[ii*dNF+jj];
    }
        cout<<"";
    }
    cout<<"   ";
}

template <class T> void quad<T>::print_local_sd()
//void quad::print_Eval_local_sd()
{
    cout<<" Shape Funtion Derivatives values: ";
    for (auto ii=dsF_vec.begin();ii!=dsF_vec.end();++ii){
        cout<<" "<< *ii;
    }
    cout<<" Shape Funtion Derivatives values: "<<local_sd_vec.size()<<" "<<local_sd_vec[0].size();
    for (auto ii=0;ii<local_sd_vec.size();++ii){
        cout<<"";
    for (auto jj=0;jj<local_sd_vec[0].size();++jj){
        cout<<" "<< local_sd_vec[ii][jj];
    }
    }
    cout<<" Shape Funtion Eta Derivatives values: ";
    for (auto ii=dsFdeta_vec.begin();ii!=dsFdeta_vec.end();++ii){
        cout<<" "<< *ii;
    }
    cout<<" Shape Funtion Xi Derivatives values: ";
    for (auto ii=dsFdxi_vec.begin();ii!=dsFdxi_vec.end();++ii){
        cout<<" "<< *ii;
    }
    cout<<"   ";
}

template <class T> void quad<T>::Eval_Jm2D(T x[NN], T y[NN]){
    int ii,jj;
    //
    // Initialization
    for (ii=0;ii<ND;++ii){
    for (jj=0;jj<ND;++jj){
        Jm2D[ii][jj] = 0.0;
    }            
    }    
    //
    if(Jm2D_vec.empty()!=true){
    while (Jm2D_vec.size() > 0) {
        Jm2D_vec.pop_back();
    }
    }
    //
    //for (auto ii=dsFdeta_vec.begin();ii!=dsFdeta_vec.end();++ii){
    for (ii=0;ii<NF;++ii){
        Jm2D[0][0] +=  dsF_vec[ii]*x[ii];
        Jm2D[1][0] +=  dsF_vec[ii]*y[ii];
        jj = ii+NF;
        Jm2D[0][1] += dsF_vec[jj]*x[ii];
        Jm2D[1][1] += dsF_vec[jj]*y[ii];
        //coord_vec.push_back(x[ii]);
        //coord_vec.push_back(y[ii]);
    }
    for(int ii=0;ii<ND;ii++){
    for(int jj=0;jj<ND;jj++){
        Jm2D_vec.push_back(Jm2D[ii][jj]);
    }        
    }
}


template <class T> void quad<T>::Eval_Jm2D(vector<T> xy){
    int ii,jj,iex,iey;
    //
    // Initialization
    for (ii=0;ii<ND;++ii){
    for (jj=0;jj<ND;++jj){
        Jm2D[ii][jj] = 0.0;
    }            
    }    
    //
    if(Jm2D_vec.empty()!=true){
    while (Jm2D_vec.size() > 0) {
        Jm2D_vec.pop_back();
    }
    }
    //
    //for (auto ii=dsFdeta_vec.begin();ii!=dsFdeta_vec.end();++ii){
    for (ii=0;ii<NF;++ii){
        iex = ND*(ii);
        iey = ND*(ii)+1;
        Jm2D[0][0] +=  dsF_vec[ii]*xy[iex];
        Jm2D[1][0] +=  dsF_vec[ii]*xy[iey];
        jj = ii+NF;
        Jm2D[0][1] += dsF_vec[jj]*xy[iex];
        Jm2D[1][1] += dsF_vec[jj]*xy[iey];
        //coord_vec.push_back(x[ii]);
        //coord_vec.push_back(y[ii]);
    }
    for(int ii=0;ii<ND;ii++){
    for(int jj=0;jj<ND;jj++){
        Jm2D_vec.push_back(Jm2D[ii][jj]);
    }        
    }
}

template <class T> void quad<T>::print_Jm2D(){
//void quad::print_J(){
    cout<<" Jacobian Matrix Jm2D: ";
    for (int ii=0;ii<ND;ii++){
    for (int jj=0;jj<ND;jj++){
        cout<<" "<< Jm2D[ii][jj];
    }
    cout<<""<<endl;
    }
    cout<<" Jacobian Matrix Jm2D_vec: ";
    for (int ii=0;ii<ND;ii++){
    for (int jj=0;jj<ND;jj++){
        cout<<" "<< Jm2D_vec[ii*ND+jj];
    }
    cout<<""<<endl;
    }
    cout<<"   ";
}

template <class T> void quad<T>::print_invJm2D(){
    cout<<" Inverse of Jacobian Matrix: ";
    for (int ii=0;ii<ND;ii++){
    for (int jj=0;jj<ND;jj++){
        cout<<" "<< invJm2D[ii][jj];
    }
    cout<<""<<endl;
    }
    cout<<"   ";
}

// Function for calculating the global Cartesian derivatives of the shape functions
template <class T> void quad<T>::Eval_global_sd()
{
    if(J_flag == 0)
        cout<<" WARNINN! : JacobianMatrix is singular";
    //for(int jj = 0; jj<invJm2D[0].size();jj++){
    for(int ii = 0; ii<ND;ii++){
        for(int jj = 0; jj<NF;jj++){
             //temp[ii][jj] = invJm2D[ii,jj]*local_sd_vec[ii,jj];
            for(int kk = 0; kk<ND;kk++){
             //cout<<""<<invJm2D[ii][kk]<<"*"<<local_sd_vec[kk][jj]<<"+";
             global_sd[ii][jj] += invJm2D[ii][kk] * local_sd_vec[kk][jj];
        }
        ////cout<<"= "<<global_sd[ii][jj]<<"";
        //cout<<" "<<global_sd[ii][jj];
        global_sd_vec.push_back(global_sd[ii][jj]);
        }
        ////cout<<" ========";
        //cout<<"";
    }
}

// Function for printing the global Cartesian derivatives of the shape functions
template <class T> void quad<T>::print_global_sd()
{
    if(J_flag == 0)
        cout<<" WARNINN! : JacobianMatrix is singular";
    for(int ii = 0; ii<ND;ii++){
        for(int jj = 0; jj<NF;jj++){
        ////cout<<"= "<<global_sd[ii][jj]<<"";
        cout<<" "<<global_sd[ii][jj];
        }
        ////cout<<" ========";
        cout<<"";
    }
}

// Function to evaluate the B matrix in 2D neccesary for computing strains
template <class T> void quad<T>::Eval_Bm2D()
{
    int ii = 0;
    for(int jj = 0; jj<dNF;jj++){
    ////cout<<"= "<<global_sd[ii][jj]<<"";
    (jj%2!=0) ? bm2D[ii][jj] = 0.0 : bm2D[ii][jj] = global_sd[ii][int(jj/2)] ;
    // jj = 0 -> bm2D[ii][jj] = global_sd[ii][int(jj/2)]
    //jj = 1 -> bm2D[ii][jj] = 0
    (jj%2!=0) ? B_vec.push_back(0.0) : B_vec.push_back(global_sd[ii][int(jj/2)]) ;
    }
    //
    ii = 1;
    for(int jj = 0; jj<dNF;jj++){
    ////cout<<"= "<<global_sd[ii][jj]<<"";
    (jj%2!=0) ?  bm2D[ii][jj] = global_sd[ii][int(jj/2)] : bm2D[ii][jj] = 0.0;
    // jj = 0 -> bm2D[ii][jj] = 0
    //jj = 1 -> bm2D[ii][jj] = global_sd[ii][int(jj/2)]
    (jj%2!=0) ?  B_vec.push_back(global_sd[ii][int(jj/2)]) : B_vec.push_back(0.0);
    }
    //
    ii = 2;
    for(int jj = 0; jj<dNF;jj++){
    ////cout<<"= "<<global_sd[ii][jj]<<"";
    (jj%2!=0) ?  bm2D[ii][jj] = global_sd[0][int(jj/2)] : bm2D[ii][jj] = global_sd[1][int(jj/2)] ;
    // jj = 0 -> bm2D[ii][jj] = global_sd[1][int(jj/2)]
    //jj = 1 -> bm2D[ii][jj] = global_sd[0][int(jj/2)]
    (jj%2!=0) ?  B_vec.push_back( global_sd[0][int(jj/2)] ) : B_vec.push_back( global_sd[1][int(jj/2)] );
    }
    /*
    (jj%2!=0) ? bm2D[ii][jj] = 0.0 : bm2D[ii][jj] = global_sd[ii][int(jj/2)] ;
    // jj = 0 -> bm2D[ii][jj] = global_sd[ii][int(jj/2)]
    //jj = 1 -> bm2D[ii][jj] = 0
    (jj%2!=0) ? B_vec.push_back(0.0) : B_vec.push_back(global_sd[ii][int(jj/2)]) ;
    }
    //
    ii = 1;
    for(int jj = 0; jj<dNF;jj++){
    ////cout<<"= "<<global_sd[ii][jj]<<"";
    (jj%2!=0) ?  bm2D[ii][jj] = global_sd[ii][int(jj/2)] : bm2D[ii][jj] = 0.0;
    // jj = 0 -> bm2D[ii][jj] = 0
    //jj = 1 -> bm2D[ii][jj] = global_sd[ii][int(jj/2)]
    (jj%2!=0) ?  B_vec.push_back(global_sd[ii][int(jj/2)]) : B_vec.push_back(0.0);
    }
    //
    for(int jj = 0; jj<dNF;jj++){
    ////cout<<"= "<<global_sd[ii][jj]<<"";
    (jj%2!=0) ?  bm2D[ii][jj] = global_sd[0][int(jj/2)] : bm2D[ii][jj] = global_sd[1][int(jj/2)] ;
    // jj = 0 -> bm2D[ii][jj] = global_sd[1][int(jj/2)]
    //jj = 1 -> bm2D[ii][jj] = global_sd[0][int(jj/2)]
    (jj%2!=0) ?  B_vec.push_back( global_sd[0][int(jj/2)] ) : B_vec.push_back( global_sd[1][int(jj/2)] );
    }
    */
}

// Function for printing the global Cartesian derivatives of the shape functions
template <class T> void quad<T>::print_bm2D()
{
    cout<<"";
    for(int ii = 0; ii<NS;ii++){
        for(int jj = 0; jj<dNF;jj++){
        cout<<" "<<bm2D[ii][jj];
        }
        cout<<"";
    }
    //
    cout<<"";
    for(int ii = 0; ii<NS;ii++){
        for(int jj = 0; jj<dNF;jj++){
        cout<<" "<<B_vec[ii*dNF+jj];
        }
        cout<<"";
    }
}

// Function for evaluating matrix D
template <class T> void quad<T>::Eval_D_vec(T E, T nu)
{
T coef;
coef = E/(1-nu*nu);
D_vec.push_back(coef);          //11
D_vec.push_back(coef*nu);       //12
D_vec.push_back(0);             //13
D_vec.push_back(coef*nu);       //21
D_vec.push_back(coef);          //22
D_vec.push_back(0);             //23
D_vec.push_back(0);             //31
D_vec.push_back(0);             //32
D_vec.push_back(coef*(1-nu)/2); //33
}


// computes principal stresses and angle. 
template <class T> void quad<T>::Eval_principalStress(qbMatrix2<T> estress){
    //estress = [sxx, syy, tauxy]
    // pr1: principal stress along the new x axis. units are [N/m^2]
    // pr2: principal stress along the new y axis
    // alphas units are radians
    T tol = 1.0e-9; // tolerance
    T sxx,syy,sxy, fact, mean,zcheck;
    sxx = estress.GetElement(0,0);
    syy = estress.GetElement(1,0);
    sxy = estress.GetElement(2,0);
    //
    fact = sqrt( pow(0.5*(sxx - syy),2) + pow(sxy,2) );
    mean = 0.5*( sxx + syy );
    this->pr1 = mean + fact;
    this->pr2 = mean - fact;
    zcheck = fabs(this->pr1 - this->pr2); // checking if the denominator is close to zero
    if (zcheck < tol){
        // if the principal stresses are close to being equal that means an angle f pure shear
        this->alpha = 3.1416/4;
    }
    else{
        this->alpha = 0.5*atan( 2*sxy / ( this->pr1 - this->pr2 ) );
    }
    //
    //
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//https://www.geeksforgeeks.org/adjoint-inverse-matrix/
// Function to get cofactor of A[p][q] in temp[][]. n is
// current dimension of A[][]
template <class T> void quad<T>:: getCofactor(T A[ND][ND], T temp[ND][ND], int p, int q,int n)
//void quad:: getCofactor(int A[ND][ND], int temp[ND][ND], int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++) {
        for (int col = 0; col < n; col++) {
            //  Copying into temporary matrix only those
            //  element which are not in given row and
            //  column
            if (row != p && col != q) {
                temp[i][j++] = A[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1) {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
 
/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
template <class T> T quad<T>::determinant(T A[ND][ND], T n)
//int quad:: determinant(int A[ND][ND], int n)
{
    T D = 0; // Initialize result
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];
 
    T temp[ND][ND]; // To store cofactors
 
    T sign = 1; // To store sign multiplier
 
    // Iterate for each element of first row
    for (int f = 0; f < n; f++) {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint of A[ND][ND] in adj[ND][ND].
template <class T> void quad<T>:: adjoint(T A[ND][ND], T adj[ND][ND])
//void quad:: adjoint(int A[ND][ND], int adj[ND][ND])
{
    if (ND == 1) {
        adj[0][0] = 1;
        return;
    }
 
    // temp is used to store cofactors of A[][]
    T sign = 1; 
    T temp[ND][ND];
 
    for (int i = 0; i < ND; i++) {
        for (int j = 0; j < ND; j++) {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, ND);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i + j) % 2 == 0) ? 1 : -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign) * (determinant(temp, ND - 1));
        }
    }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
template <class T> bool quad<T>::inverse(T A[ND][ND], T inverse[ND][ND])
//bool quad:: inverse(int A[ND][ND], float inverse[ND][ND])
{
    // Find determinant of A[][]
    T det = determinant(A, ND);
    if (det == 0) {
        cout << "Singular matrix, can't find its inverse";
        return false;
    }
 
    // Find adjoint
    T adj[ND][ND];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) =
    // adj(A)/det(A)"
    for (int i = 0; i < ND; i++)
        for (int j = 0; j < ND; j++)
            inverse[i][j] = adj[i][j] / float(det);
 
    return true;
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//deconstructor
template <class T> quad<T>::quad::~quad()
{

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//constructor
template <class T> quadModel<T>::quadModel()
{


}

template <class T> quadModel<T>::quadModel(string elm_data_FileName)
{
string dummyLine; // auxiliary variable for catching heathers
string dummyWord; // auxiliary variable for catching words
T E_val;
T nu_val;
T mu_val;
T gt_val;
int selectFlag_val;
//
long int ne_val = 0; // number of elements
long int nds_val = 0; // number of nodes
long int nf_val = 0; // number of fixed nodes
long int ndsd_val = 0; // number of dynamic constrains
int modelType_val;
int elementType; // type of element, 4: quad linear
long int nds1ID, nds2ID,nds3ID, nds4ID, fg;
T xCoord, yCoord; // temporal node coordinates
T uval; // temporal node displacements
//
long int nodeID;
long int elementID;
//
cout << "Reading Data File";
ifstream elm_data(elm_data_FileName);
if(elm_data.is_open()){
    getline(elm_data,dummyLine);
    cout<<"File info: " << dummyLine << endl;
    elm_data >> dummyWord >> dummyWord >> E_val >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> nu_val >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> mu_val >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> gt_val >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> selectFlag_val >> dummyWord; 
    this->E = E_val;
    this->nu = nu_val;
    this->mu = mu_val;
    this->gt = gt_val;
    this->selectFlag = selectFlag_val;
    //
    //
    elm_data >> dummyWord >> dummyWord >> ne_val >> dummyWord;// dummyWord >> dummyWord >> dummyWord >> dummyWord >> ne >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> nds_val >> dummyWord;// dummyWord >> dummyWord >> dummyWord >> dummyWord >> nds >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> nf_val >> dummyWord;// dummyWord >> dummyWord >> dummyWord >> dummyWord >> dummyWord >> nf >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> ndsd_val >> dummyWord;// dummyWord >> dummyWord >> dummyWord >> dummyWord >> dummyWord >> ndsd >> dummyWord; 
    elm_data >> dummyWord >> dummyWord >> modelType_val >> dummyWord;
    this->ne = ne_val;
    this->nds = nds_val;
    this->nf = nf_val;
    this->ndsd = ndsd_val;
    this->modelType = modelType_val;
    this->nev = this->nds*ND; // total number of degrees of freedom in x and y axis
    this->nevf = this->nds*NS; // total number of degrees of freedom in x, y and xy axis
    //
    cout<<"Mechanical properties: " << endl;
	cout << "Young M." << setw(10) << "nu" << setw(15) << "mu" << setw(20)<<"h" << setw(20)<< "selectFlag" << endl;
	cout << E << setw(15) << nu << setw(15) << mu << setw(18)<<gt_val << setw(15) << selectFlag << endl;
    //
	cout << "Number Of Elements" << setw(20) << "Number Of Nodes" << endl;
	cout << this->ne << setw(25) << this->nds <<endl;
    cout << "Number Of Fixed Nodes" << setw(40)<< "Number Of Dynamic Constrains" << endl;
    cout << this->nf << setw(40)<< this->ndsd <<  endl;
    cout << "Model Type" << setw(40)<< "Global Element Thickness" << endl;
    cout << this->modelType << setw(40)<< gt <<  endl;
    cout << "Total number of degrees of freedom" << endl;
    cout << this->nev << endl;
    //
    elm_data >> dummyWord;
    cout<<"Reading element map: " << endl;
    cout<<dummyWord<< endl;
    elementMap.Resize(ne,modelType+1);
    elementMap.SetToZero();
    for(long int elementID = 0; elementID < this->ne; elementID++){
    elm_data >> elementType >> dummyWord>> nds1ID >> dummyWord>> nds2ID >> dummyWord>> nds3ID >> dummyWord>> nds4ID >> dummyWord;
    elementMap.SetElement(elementID,0,elementType);
    elementMap.SetElement(elementID,1,nds1ID-1);
    elementMap.SetElement(elementID,2,nds2ID-1);
    elementMap.SetElement(elementID,3,nds3ID-1);
    elementMap.SetElement(elementID,4,nds4ID-1);
    }
    elementMap.PrintMatrix();     
    //
    elm_data >> dummyWord;
    cout<<"Reading node coordinates map: " << endl;
    cout<<dummyWord<< endl;
    gXYvec.Resize(this->nds,ND);
    gXYvec.SetToZero();
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
    elm_data >> xCoord>> dummyWord >> yCoord >> dummyWord;
    gXYvec.SetElement(nodeID,0,xCoord);
    gXYvec.SetElement(nodeID,1,yCoord);
    }
    gXYvec.PrintMatrix();  
    //
    //
    //
    switch (selectFlag){
    //
    case 0: 
        elm_data >> dummyWord;
        cout<<"Reading node displacement map: " << endl;
        cout<<dummyWord<< endl;
        gu.Resize(this->nev,1);
        gu.SetToZero();
        for(nodeID = 0; nodeID < this->nev; nodeID++){
        elm_data >> uval >> dummyWord;
        gu.SetElement(nodeID,0,uval);
        }
        gu.PrintMatrix();     
        break;
    //
    case 1:
        elm_data >> dummyWord;
        cout<<"Reading global nodal indexes map where constrains are imposed: " << endl;
        nofElementMap.Resize(this->nf,1);
        nofElementMap.SetToZero();
        for(elementID = 0; elementID < this->nf; elementID++){
        elm_data >> nds1ID  >> dummyWord;
        nofElementMap.SetElement(elementID,0,nds1ID-1);// nodes enumeration goes fron 0 to nds
        }
        nofElementMap.PrintMatrix();             
        //
        elm_data >> dummyWord;
        cout<<"Reading global nodal map of constrains values: 1: fixed node, 0: free node " << endl;
        ifp.Resize(this->nf,2);
        ifp.SetToZero();
        for(nodeID = 0; nodeID < this->nf; nodeID++){
        elm_data >> xCoord >> dummyWord >> yCoord >> dummyWord;
        ifp.SetElement(nodeID,0,xCoord);
        ifp.SetElement(nodeID,1,yCoord);
        }
        ifp.PrintMatrix();
        //
        //
        // prescribed displacements
        pre.Resize(this->nf,2);
        pre.SetToZero();
        //
        //
        elm_data >> dummyWord;
        cout<<"Reading global nodal indexes map where forces are imposed: " << endl;
        fxyElementMap.Resize(this->ndsd,1);
        fxyElementMap.SetToZero();
        for(elementID = 0; elementID < this->ndsd; elementID++){
        elm_data >> nds1ID  >> dummyWord;
        fxyElementMap.SetElement(elementID,0,nds1ID-1); // nodes enumeration goes fron 0 to nds
        }
        fxyElementMap.PrintMatrix();             
        //
        elm_data >> dummyWord;
        cout<<"Reading global nodal map of force values " << endl;
        //cout<<dummyWord<< endl;
        //fxy.Resize(ndsd,2);
        fxy.Resize(this->nev,1);
        fxy.SetToZero();
        for(nodeID = 0; nodeID < this->ndsd; nodeID++){
        elm_data >> xCoord >> dummyWord >> yCoord >> dummyWord;
        fg = fxyElementMap.GetElement(nodeID,0);
        //cout << fg << " " << fg*ND+0 << " "<< fg*ND+1<< endl;
        fxy.SetElement(fg*ND+0,0,xCoord);
        fxy.SetElement(fg*ND+1,0,yCoord);
        }
        fxy.PrintMatrix(); 
        //
        //     
        break;
    //
    default:
        cout<< "No Selection"<< endl;
        break;
    }

}
else{
    cout<< "Error! Input file: "<<elm_data_FileName<<"could not be opened.WARNING! USING DEAFAULT VALUES";
    /*
    //warning: overflow in implicit constant conversion
    E = 2e11;
    nu = 0.25;
    mu = 8e10;
    h = 10;
    selectFlag = 0;
    */
}


// Initializing nodal coordinates
eXYvec1.push_back(0.0);
eXYvec1.push_back(0.0);
eXYvec1.push_back(100.0);
eXYvec1.push_back(0.0);
eXYvec1.push_back(100.0);
eXYvec1.push_back(100.0);
eXYvec1.push_back(0.0);
eXYvec1.push_back(100.0);

//
eXYvec2.push_back(0.0);
eXYvec2.push_back(0.0);
eXYvec2.push_back(100.0);
eXYvec2.push_back(0.0);
eXYvec2.push_back(100.0);
eXYvec2.push_back(100.0);
eXYvec2.push_back(0.0);
eXYvec2.push_back(100.0);




//////////////////////
//
// Initializing displacements

//
ue_xy1.push_back(0.000); //ux1
ue_xy1.push_back(0.000); //uy1
ue_xy1.push_back(0.045); //ux2
ue_xy1.push_back(0.004); //uy2
ue_xy1.push_back(0.046); //ux3
ue_xy1.push_back(0.027); //uy3
ue_xy1.push_back(0.008); //ux4
ue_xy1.push_back(0.029); //uy4
//
// Uniaxial stress
ue_xy2.push_back(0.000); //ux1
ue_xy2.push_back(0.000); //uy1
ue_xy2.push_back(0.100); //ux2
ue_xy2.push_back(0.000); //uy2
ue_xy2.push_back(0.100); //ux3
ue_xy2.push_back(-0.03); //uy3
ue_xy2.push_back(0.000); //ux4
ue_xy2.push_back(-0.03); //uy4


}

template <class T> quadModel<T>::quadModel(T E_val,T nu_val, T t_val)
{
//cout << "Creating Element";
this->E = E_val;
this->nu = nu_val;
this->gt = t_val;
}

template <class T> vector<T> quadModel<T>::get_sF_vals(T xi, T eta)  //returns local shape function values
{
    quad<T> quad1(xi,eta);
    quad1.Eval_local_sf(xi,eta);
    //quad1.print_sf();
    return quad1.get_local_sf_vec();
}

template <class T> qbMatrix2<T> quadModel<T>::get_dsF_matrix(T xi, T eta)  //returns local shape function values
{
    quad<T> quad1(xi,eta);
    quad1.Eval_local_sd(xi,eta);
    //quad1.print_sf();
    qbMatrix2<T> temp(ND,NF,quad1.get_dsF_vec());
    return temp;
}

 //basic modeller to compute stress and strain in a quad element
template <class T> void quadModel<T>::basicAssembler()
{
    T xi, eta;
    T J1Det;
    //qbMatrix2<T> Bm2D(NS,NEV);
    //
    //
    xi  = 0.25;
    eta = 0.50;
    //
    // Creating quad element
    quad<T> quad1(xi,eta);
    //
    cout << "General Plane Stres-Strain case: ";
    // Computing b matrix
    qbMatrix2<T> Bm2D1 = this->Bmtx_quad(xi,eta,eXYvec1,J1Det);
    //qbMatrix2<T> Bm2D = Bmtx_quad(xi,eta);
    //Bm2D.PrintMatrix();
    //
    // Computing D matrix
    quad1.Eval_D_vec(E, nu);
    qbMatrix2<T> Dm2D1(NS,NS,quad1.get_D_vec());
    //D.PrintMatrix();
    //
    // Computing Strain
    cout << "Printing local Strain: ";
    qbMatrix2<T> ue1(NEV,1,ue_xy1);
    qbMatrix2<T> stn1 = Bm2D1*ue1;
    stn1.PrintMatrix();
    //
    // Computing Stress
    cout << "Printing local Stress: ";
    qbMatrix2<T> str1 = Dm2D1*stn1;
    str1.PrintMatrix();    
    //
    //
   cout << "General Uniaxial Stress case: ";
    // Computing b matrix
    qbMatrix2<T> Bm2D2 = this->Bmtx_quad(xi,eta,eXYvec2,J1Det);
    //qbMatrix2<T> Bm2D = Bmtx_quad(xi,eta);
    //Bm2D.PrintMatrix();
    //
    // Computing D matrix
    quad1.Eval_D_vec(E, nu);
    qbMatrix2<T> Dm2D2(NS,NS,quad1.get_D_vec());
    //D.PrintMatrix();
    //
    // Computing Strain
    cout << "Printing local Strain: ";
    qbMatrix2<T> ue2(NEV,1,ue_xy2);
    qbMatrix2<T> stn2 = Bm2D2*ue2;
    stn2.PrintMatrix();
    //
    // Computing Stress
    cout << "Printing local Stress: ";
    qbMatrix2<T> str2 = Dm2D2*stn2;
    str2.PrintMatrix();    
}

 //basic modeller to compute test local stiffness matrix
template <class T> void quadModel<T>::basicElementAssembler()
{
    cout<<"Basic Element Assembler "<< endl;
    long int nds1; // temporal node id
    long int elementID,nodeID,nodesPerElement,id;
    // elementID: element index
    // nodeID: node index
    // nodesPerElement : number of nodes in element "elementID"
    //
    // Computing D matrix
    this->Eval_D_vec(this->E,this->nu);
    qbMatrix2<T> Dm2D1(NS,NS,this->get_D_vec());
    //cout << "Printing D matrix Stress: ";
    //Dm2D1.PrintMatrix(); 
    //
    // Sweeping elements
    for(elementID = 0; elementID < this->ne; elementID++){
    // seeping nodes for Initializing nodal coordinates
        nodeID = 0;
        nodesPerElement = this->elementMap.GetElement(elementID,nodeID);
        vector<T> exyVec1; // local vector of nodal corrdinates
        vector<T> euVec1; // local vector of nodal displancements
        //
        for(nodeID = 1; nodeID < nodesPerElement+1; nodeID++){
            nds1 = this->elementMap.GetElement(elementID,nodeID);
            for(id = 0; id < ND; id++){
                exyVec1.push_back(this->gXYvec.GetElement(nds1,id)); // id-th coordinate
                euVec1.push_back(this->gu.GetElement(ND*(nds1)+id,0)); // displancement on the id-th axis
            }
            //cout << nodeID << " " <<nds1 << endl;
        }
        //
        //qbMatrix2<T> XY1(NF,ND,xye1);
        //XY1.PrintMatrix();
        qbMatrix2<T> eu1(NEV,1,euVec1);
        //eu1.PrintMatrix();
        cout << "Stiff Matrix for element "<< elementID << endl;
        qbMatrix2<T> Ke = this -> Ematx_quad(NEV,NG,exyVec1, Dm2D1,this->gt);
        Ke.PrintMatrix(); 
        cout << "Forces acting on element "<< elementID << endl;    
        qbMatrix2<T> Fe = Ke*eu1;
        Fe.PrintMatrix();
    }
    //
//
}


 //basic modeller to compute global Kg stiffnes matrix
template <class T> void quadModel<T>::basicKgAssembler()
{
    cout<<"Basic Kg Assembler "<< endl;
    long int nds1; // temporal node id
    long int elementID,nodeID,nodesPerElement;
    long int ie,ig,je,jg,id,jd;
    int iin,jjn;
    // elementID: element index
    // nodeID: node index
    // nodesPerElement : number of nodes in element "elementID"
    //
    Kg.Resize(nev,nev);
    //Kg.SetToZero();
    //Kg.PrintMatrix();
    //
    //
    // Computing D matrix
    this->Eval_D_vec(this->E,this->nu);
    qbMatrix2<T> Dm2D1(NS,NS,this->get_D_vec());
    //
    // Sweeping elements
    for(elementID = 0; elementID < this->ne; elementID++){
    // sweeping nodes for Initializing nodal coordinates
        nodeID = 0;
        nodesPerElement = this->elementMap.GetElement(elementID,nodeID);
        vector<T> exyVec1; // local vector of nodal corrdinates
        //
        for(nodeID = 1; nodeID < nodesPerElement+1; nodeID++){
            nds1 = this->elementMap.GetElement(elementID,nodeID);
            for(id = 0; id < ND; id++){
            exyVec1.push_back(this->gXYvec.GetElement(nds1,id)); // id-th coordinate
            }
        }
        //
        //cout << "Stiff Matrix for element "<< elementID << endl;
        qbMatrix2<T> Ke = this -> Ematx_quad(NEV,NG,exyVec1, Dm2D1,this->gt);
        //Ke.PrintMatrix(); 
        //
        //
        // Transfering elements from local stiffness matrix Ke to global stiffness matrix Kg
        // sweeping i nodes
        for(iin = 1; iin < nodesPerElement+1; iin++){
            for(id = 0; id < ND; id++){
                ie = (iin-1)*ND + id;
                ig = (this->elementMap.GetElement(elementID,iin))*ND + id;
        // sweeping j nodes
        for(jjn = 1; jjn < nodesPerElement+1; jjn++){
            for(int jd = 0; jd < ND; jd++){
                je = (jjn-1)*ND + jd;
                jg = (this->elementMap.GetElement(elementID,jjn))*ND + jd;
                T tempVal = Kg.GetElement(ig,jg) + Ke.GetElement(ie,je);
                Kg.SetElement(ig,jg,tempVal);
            }
        }        // sweeping j nodes
            }
        }        // sweeping i nodes
        
    }// Sweeping elements
    //
    //cout<<" Printing Kg Matrix:"<<endl;
    //Kg.PrintMatrix();
    //
//
}



 //basic modeller to compute global RHS vector
template <class T> void quadModel<T>::basicRhsAssembler()
{
    //
    long int iv,id,ii,ig;
    T tempVal = 0.0;
    //
    cout<<"Basic RHS Assembler "<< endl;
    //
    // lets copy current RHS vector and Kg stiffness matrices
    //qbMatrix2<T> Kgm  = Kg;
    //qbMatrix2<T> RHSm = fxy;
    Kgm  = Kg;
    RHSm = fxy;
    //RHSm.SetElement(10,0,0);
    //cout<<"Original RHS "<< endl;
    //fxy.PrintMatrix();
    //cout<<"Back Up RHS "<< endl;
    //RHSm.PrintMatrix();
    //cout<<"Original RHS "<< endl;
    //fxy.PrintMatrix();
    //
    //
    /* lets modify the RHSM for fixed degrees of freedom (dof), then lets set to zero the associated rows and 
        columns in the stiffness matrix Kgm.*/
    //
    // sweeping fixed nodes
    for(iv = 0; iv < nf; iv++){
        // sweeping node degrees of fredom
        for(id = 0; id < ND; id++){
            // checking if the node is fixed
            if(ifp.GetElement(iv,id)== int(1) ){
                //
                ig = nofElementMap.GetElement(iv,0)*ND+id;
                //cout << iv << " " << nofElementMap.GetElement(iv,0) << " " <<  ig << " "<< id << " " << ifp.GetElement(iv,id) << endl;               
                for(ii = 0; ii < nev; ii++){ // sweeping degrees of freedom
                    tempVal = RHSm.GetElement(ii,0) -  pre.GetElement(iv,id)*Kg.GetElement(ii,ig);
                    RHSm.SetElement(ii,0,tempVal);
                    Kgm.SetElement(ig,ii, 0.0);
                    Kgm.SetElement(ii,ig, 0.0);                    
                } // sweeping degrees of freedom
                Kgm.SetElement( ig,ig,Kg.GetElement(ig,ig) );
                //
            }   // checking if the node is fixed
        }   // sweeping node degrees of fredom
    }    // sweeping fixed nodes
    //
    //
    //After that lets set the diagoan term in Kgm associate with that dof to its original value in Kg 
    // sweeping fixed nodes    
    for(iv = 0; iv < nf; iv++){
        // sweeping node degrees of fredom
        for(id = 0; id < ND; id++){
            // checking if the node is fixed
            if(ifp.GetElement(iv,id)== int(1) ){
                //
                ig = nofElementMap.GetElement(iv,0)*ND+id;
                tempVal = pre.GetElement(iv,id)*Kg.GetElement(ig,ig);
                RHSm.SetElement(ig,0,tempVal);
                //
            }   // checking if the node is fixed
        }   // sweeping node degrees of fredom
    }    // sweeping fixed nodes    
    //
    //
//
//
}


template <class T> void quadModel<T>::basicKgmRHSmSolver()     //returns the global displacements ug after solving the system [Kgm]*ug = RHSm   
{
    /*
    this method inverts the matrix using the LU Factoriazation defined in class qbmatrix
    part of the qbLinAlg linear algebra library
    Copyright (c) 2021 Michael Bennett
    MIT license
    */
//cout<<"Kg stiffness matrix"<<endl;
//Kgm.PrintMatrix();
//cout<<"RHSm"<<endl;;
//RHSm.PrintMatrix();
//qbMatrix2<T> u1(nev,1);
//u1.SetToZero();
//
// Computing inverse matrix
qbMatrix2<T> Kg_inv = Kgm; // initializing inverse matrix
if(Kg_inv.Inverse()){
    //cout<<"Invers Matrix"<<endl;
    //Kg_inv.PrintMatrix(10)<<endl;
    //u1 = Kg_inv*RHSm;
    cout<<"Global Displacements "<<endl;
    //u1.PrintMatrix();
    gu.Resize(this->nev,1);
    gu.SetToZero();
    gu = Kg_inv*RHSm;
    cout << endl;
    gu.PrintMatrix();
}
}

template <class T> void quadModel<T>::qbLinSolveKgmRHSmSolver()     //returns the global displacements ug after solving the system [Kgm]*ug = RHSm   
{
    /*
    this method inverts the matrix using the Gauss-Seidel substitution,
    part of the qbLinAlg linear algebra library
    Copyright (c) 2021 Michael Bennett
    MIT license
    */
   //
   int infoFlag; // indicates if the solver maneg to find a solution
   //
    //cout<<"Kg stiffness matrix"<<endl;
    //Kgm.PrintMatrix()<<endl;
    //cout<<"RHSm"<<endl;
    qbVector<T> RHSmVec(this->nev); // global nodal map of force values
    //RHSmVec = RHSm.vector();
    for(int ii=0;ii<nev;ii++)
        RHSmVec.SetElement(ii,RHSm.GetElement(ii,0));
    //RHSmVec.PrintVector()<<endl;
    //
    qbVector<T> uvec1(this->nev);
    infoFlag = qbLinSolve(Kgm, RHSmVec, uvec1);
    cout<<"Global Displacements " << infoFlag << endl;
    //uvec1.PrintVector();
    gu.Resize(this->nev,1);
    gu.SetToZero();
    for(int ii = 0; ii < this->nev; ii++){
    gu.SetElement(ii,0,uvec1.GetElement(ii));
    }
    cout << endl;
    gu.PrintMatrix();   
}


    // dgeev_ is a symbol in the LAPACK library files
    //extern "C" {
    //    //dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    //            //dgesv(N,      NRHS,   A,          LDA,    IPIV,   B,          LDB,    INFO)
    //extern int dgesv_(long int*,int*,long double*,long int*,int*,long double*,long int*, int*);
    //}


template <class T> void quadModel<T>::generalLAPCAKSolveKgmRHSmSolver()     //returns the global displacements ug after solving the system [Kgm]*ug = RHSm   
{
    /*
    this method solves the system with a general LAPACK solver
    */
   //
   long int N = this->nev, LDA = N, LDB = N;
   //long int N = 3, LDA = N, LDB = N;
   //int N = 3, LDA = N, LDB = N;
   int NRHS = 1;
   int IPIV[N]; // vector for storing the pivot operations
   int INFO = 0; // indicates if the solver maneg to find a solution
   //
   //cout<<"Kg stiffness matrix"<<endl;
      //long double A[this->nev*this->nev];
      double A[N*N];/* = {
          6.80, -2.11, 5.66,
          -6.05, -3.30, 5.36,
          -0.45, 2.58, -2.70
      };*/
      for (int ii = 0; ii < N; ii++){
        for (int jj = 0; jj < N; jj++){
          A[ii*N+jj] = Kgm.GetElement(ii,jj);
          //cout << A[ii*N+jj] << " ";
        }
        //cout << endl;
      }
    //Kgm.PrintMatrix();
    //cout<<"RHSm"<<endl;
      double B[N];
      //double B[N] = {4.02, 6.19, -8.22};
      for (int ii = 0; ii < N; ii++){
          B[ii] = RHSm.GetElement(ii,0);
          //cout << B[ii] << " ";
        //cout << endl;
      }
    //
    //dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    //INFO = dgesv_(&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);
    INFO = LAPACKE_dgesv(LAPACK_COL_MAJOR,N,NRHS,A,LDA,IPIV,B,LDB);
    cout<<"Global Displacements calculated with LAPACK DGESV " << INFO << endl;
    //qbMatrix2<T> u1(N,1,B); // global nodal map of force values
    //u1.PrintMatrix();
    gu.Resize(this->nev,1);
    gu.SetToZero();
    for(int ii = 0; ii < N; ii++){
    gu.SetElement(ii,0,B[ii]);
    }
    cout << endl;
    gu.PrintMatrix();   
}

template <class T> void quadModel<T>::bandedLAPCAKSolveKgmRHSmSolver()     //returns the global displacements ug after solving the system [Kgm]*ug = RHSm   
{
    /*
    this method solves the system with a general LAPACK solver
    */
   //
    //
    /*
   N = 4, LDB = N;
   KU = 2 ;// number of upper subdiagonals
   KL = 1;// number of lower subdiagonals
   int LDAB = 2*KL+KU+1;
   int NRHS = 1;
   int IPIV[N]; // vector for storing the pivot operations
   int INFO = 0; // indicates if the solver maneg to find a solution
    double Am[N*N] =  {-0.23  , 2.54  , -3.66     ,  0,  
                       -6.98  , 2.46  , -2.73     ,  -2.73,  
                        0     , 2.56  ,  2.46     ,  4.07,  
                        0     , 0     ,  -4.78    ,  3.82
                    };
    qbMatrix2<T>AB(N,N,Am);
    cout<<"Original  stiffness matrix";      
    AB.PrintMatrix();
    double A[LDAB*N];
    cout<<"Default Kg stiffness matrix";
      for (int ii = 0; ii < LDAB; ii++){
        for (int jj = 0; jj < N; jj++){
          A[ii*N+jj] = 0.0;
          cout << A[ii*N+jj] << " ";
        }
        cout << endl;
      }
      int count = KL*N;
      int kun = 0;
      for (int ii = N-KU; ii > -1; ii--){
        cout << "0x" << ii <<" ";
        count = count + ii;
        for (int jj= ii+kun; jj < N-kun; jj++){
          //A[(KL+KU+ii-jj)*N+jj] = Kgm.GetElement(ii,jj);
          cout <<count<< "a"<< jj-ii << ""<< jj << " ";
          A[count] = AB.GetElement(jj-ii,jj);
          count++;
        }
        cout << endl;
      }

      for (int ii = 1; ii < N-KL-1; ii++){
        for (int jj= ii+kun; jj < N-kun; jj++){
          //A[(KL+KU+ii-jj)*N+jj] = Kgm.GetElement(ii,jj);
          cout <<count << "a"<< jj << ""<< jj-ii << " ";
          A[count] = AB.GetElement(jj,jj-ii);
          count++;
        }
        cout << "0x" << ii <<" ";
        count = count + ii;
        cout << endl;
      }


        cout << endl;
        cout << endl;
      for (int ii = 0; ii < LDAB; ii++){
        for (int jj = 0; jj < N; jj++){
          cout << A[ii*N+jj] << " ";
        }
        cout << endl;
      }

    double b[4] =  {4.42, 27.13, -6.14, 10.5};
    */
   //
   //1 row of elemets
   //long int N = this->nev, LDA = N, LDB = N;
   //int KU = 7 ;// number of upper subdiagonals
   //int KL = 7;// number of lower subdiagonals
   long int N = this->nev, LDA = N, LDB = N;
   int KU = 14 ;// number of upper subdiagonals
   int KL = 14;// number of lower subdiagonals
   int LDAB = 2*KL+KU+1;
   int NRHS = 1;
   int IPIV[N]; // vector for storing the pivot operations
   int INFO = 0; // indicates if the solver maneg to find a solution
   int count = KL*N;
   //long int N = 3, LDA = N, LDB = N;
   //int N = 3, LDA = N, LDB = N;
    //cout<<"Original  stiffness matrix"<<endl;  
    //Kgm.PrintMatrix();
    double A[LDAB*N];
    //cout<<"Default Kg stiffness matrix"<<endl;
      for (int ii = 0; ii < LDAB; ii++){
        for (int jj = 0; jj < N; jj++){
          A[ii*N+jj] = 0.0;
          //cout << A[ii*N+jj] << " ";
        }
        //cout << endl;
      }
      //
      for (int ii = KU; ii > -1; ii--){
        //cout << "0x" << ii <<" ";
        count = count + ii;
        for (int jj= ii; jj < N; jj++){
          //cout <<count<< "a"<< jj-ii << ""<< jj << " ";
          A[count] = Kgm.GetElement(jj-ii,jj);
          count++;
        }
        //cout << endl;
      }
      //
      for (int ii = 1; ii < KL+1; ii++){
        for (int jj= ii; jj < N; jj++){
          //A[(KL+KU+ii-jj)*N+jj] = Kgm.GetElement(ii,jj);
          //cout <<count << "a"<< jj << ""<< jj-ii << " ";
          A[count] = Kgm.GetElement(jj,jj-ii);
          count++;
        }
        //cout << "0x" << ii <<" ";
        count = count + ii;
        //cout << endl;
      }
      //
        //cout << endl;
        //cout << endl;
      for (int ii = 0; ii < LDAB; ii++){
        for (int jj = 0; jj < N; jj++){
          //cout << A[ii*N+jj] << " ";
        }
        //cout << endl;
      }
      //
    //cout<<"RHSm"<<endl;
      double B[N];
      //double B[N] = {4.02, 6.19, -8.22};
      for (int ii = 0; ii < N; ii++){
          B[ii] = RHSm.GetElement(ii,0);
          //cout << B[ii] << " ";
        //cout << endl;
      }
      //
    //dgesv(N,NRHS,A,LDA,IPIV,B,LDB,INFO)
    //INFO = dgesv_(&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);
    //INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,N,KL,KU,NRHS,A,N,IPIV,B,LDB);
    //INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,N,1,2,2,A,N,IPIV,b,1);
    //INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, N, KL,KU, NRHS, A, N, IPIV, b, 1);
    INFO = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, N, KL,KU, NRHS, A, N, IPIV, B, 1);
    cout<<"Global Displacements calculated with LAPACK DGBSV " << INFO << endl;
    //qbMatrix2<T> u1(N,1,B); // global nodal map of force values
    //u1.PrintMatrix();
    gu.Resize(this->nev,1);
    gu.SetToZero();
    for(int ii = 0; ii < N; ii++){
    gu.SetElement(ii,0,B[ii]);
    }
    cout << endl;
    gu.PrintMatrix();   
}


// writes nodal stresss, strains, deformation, initial and final location
template <class T> void quadModel<T>::vtkData(string FileName){
    //string FileName = "vtkData.vtk";
	ofstream OutputFile(FileName);
    //
    int nne = NN; // number of nodes per element
    long int N = this->nev; // total number of degrees of freedom
    long int ig;
    //
    //Save a 3-D scalar array in VTKk format
    // filename in VTK format using the "POYGONS" scheme
    //
    ////////////////////////////////// Writing node data
    OutputFile << "# vtk DataFile Version 3.0" << endl;
    OutputFile << "2D Stress-Strain Distribution" << endl;
    OutputFile << "ASCII" << endl;
    //
    OutputFile << "DATASET POLYDATA" << endl;
    OutputFile << "POINTS " << nds << " float" << endl;
    //printf( "POINTS %d float",nds << endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        // xCoord << yCoord << zCoord
        OutputFile <<  gXYvec.GetElement(nodeID,0) << " " << gXYvec.GetElement(nodeID,1) << " " << 0.0 << endl;
    }
    //
    OutputFile << "POLYGONS " << ne << " " << ne*(nne+1) << " " << endl;
    for(long int elementID = 0; elementID < this->ne; elementID++){
    //elementMap.GetElement(elementID,0) << elementMap.GetElement(elementID,1) << " " << elementMap.GetElement(elementID,2) << " " << elementMap.GetElement(elementID,3) << " " << elementMap.GetElement(elementID,4) << " " << endl;
    OutputFile << nne << " " << elementMap.GetElement(elementID,1) << " " << elementMap.GetElement(elementID,2) << " " << elementMap.GetElement(elementID,3) << " " << elementMap.GetElement(elementID,4) << " " << endl;
    }
    //
    //
    OutputFile << "CELL_DATA "<< ne << endl;   
    OutputFile << "NORMALS cell_normals float"<< endl;
    for(long int elementID = 0; elementID < this->ne; elementID++){
        OutputFile << 0.0 << " " << 0.0 << " " << 1.0 << endl;
    }
    
    OutputFile << "POINT_DATA " << nds << endl;   
    OutputFile << "SCALARS x_displacement float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gu.GetElement(nodeID*ND+0,0) << endl;
    }
    OutputFile << "SCALARS y_displacement float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gu.GetElement(nodeID*ND+1,0) << endl;
    }
    ///////////////////////////////////////////////////Strains
    OutputFile << "SCALARS x_strain float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    //gstrain.PrintMatrix();
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        ig = (nodeID*NS)+0;
        OutputFile << gstrain.GetElement(ig,0) << endl;
    }
    OutputFile << "SCALARS y_strain float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gstrain.GetElement(nodeID*NS+1,0) << endl;
    }
    OutputFile << "SCALARS xy_strain float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gstrain.GetElement(nodeID*NS+2,0) << endl;
    }
    /////////////////////////////////////////////////Stresses
    OutputFile << "SCALARS x_stress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gstress.GetElement(nodeID*NS+0,0) << endl;
    }
    OutputFile << "SCALARS y_stress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gstress.GetElement(nodeID*NS+1,0) << endl;
    }
    OutputFile << "SCALARS xy_stress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gstress.GetElement(nodeID*NS+2,0) << endl;
    }
    ///////////////////////////////////////////////// Printcipal Stresses
    OutputFile << "SCALARS x_pstress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpstress.GetElement(nodeID*NS+0,0) << endl;
    }
    OutputFile << "SCALARS y_pstress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpstress.GetElement(nodeID*NS+1,0) << endl;
    }
    OutputFile << "SCALARS xy_pstress float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpstress.GetElement(nodeID*NS+2,0) << endl;
    }
    ///////////////////////////////////////////////// Global Gauss Points
    OutputFile << "SCALARS x_pXYvec float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpXYvec.GetElement(nodeID,0) << endl;
    }
    OutputFile << "SCALARS y_pXYvec float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpXYvec.GetElement(nodeID,1) << endl;
    }
    OutputFile << "SCALARS z_pXYvec float 1"<< endl;
    OutputFile << "LOOKUP_TABLE default"<< endl;
    for(long int nodeID = 0; nodeID < this->nds; nodeID++){
        OutputFile << gpXYvec.GetElement(nodeID,2) << endl;
    }
    //
    //
    ///////////////////////////////// Writing element data
    //for(long int elementID = 0; elementID < this->ne; elementID++){
    //    OutputFile << eStrain << endl;
    //}
    
	OutputFile.close();
}


// compute strains and stresses on elements looping over gauss integration points
template <class T> void quadModel<T>::Stresses(){
//
long int elementID,nodeID,nodesPerElement,id,ie,nds1;
long int localii,globalii;
long int N = this->nev; // total number of degrees of freedom
long int Nf = this->nevf; // total number of degrees of freedom
//
// Lets initialize the gauss points coordinate arrays
vector<T> pgp_vec; // vector of Gauss Coordinates
T c_val = 0.577350269189626; // CHECK accuracy
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(0.0);    
    //
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(0.0);    
    //
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(0.0);    
    //
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(0.0);    
    int ngp = 4; // number of Gauss Points
qbMatrix2<T> pgp(NF,3,pgp_vec); // local gauss points    
gpXYvec.Resize(this->nds,3); // global gauss points
gpXYvec.SetToZero();
//
// lets initialize global strains
gstrain.Resize(Nf,1);
gstrain.SetToZero();
//
// lets initialize global stresses
gstress.Resize(Nf,1);
gstress.SetToZero();
// lets initialize global principal stresses
gpstress.Resize(Nf,1);
gpstress.SetToZero();
//
// lets initialize local displacements
qbMatrix2<T> eu(NEV,1); //nev = nn*nd, Number of variables in an element
//qbMatrix2<T> exy(NF,ND); //nev = nn*nd, Number of variables in an element
eu.SetToZero();
qbMatrix2<T> B(NS,NEV);
qbMatrix2<T> D(NS,NS,this->get_D_vec());
qbMatrix2<T> estrain(NS,1);// = B*eu2;
qbMatrix2<T> estress(NS,1);// = D*estrain;
qbMatrix2<T> pstress(NS,1);// principal stresses and algle of the new coordinate system
qbMatrix2<T> pXYvec(1,NS);
//
// lets initialize global Gauss point counter
int lgp = 0;
//
//lets loop over all elements, ne: total number of ellements
cout << "Element Stresses" << endl;
cout << "Emt ID"  << " " << "Node ID" << " " << "Gauss Point" << " " << "GP X-coord." << " " << "GP Y-coord." << " " << "XXstress" << " " << "YYstress" << " " << "XYstress" << " " << "Pricipal Stress 1" << " " << "Pricipal Stress 2" << endl;
for(elementID = 0;elementID < this->ne; elementID++){
    //////////////////////////// lets extract element nodal coordinates and displancements
    nodeID = 0;
    nodesPerElement = this->elementMap.GetElement(elementID,nodeID);
    vector<T> exyVec1; // local vector of nodal corrdinates
    vector<T> euVec1; // local vector of nodal displancements
    //for(nodeID = 0; nodeID < NN; nodeID++){     //NN: nodes per element
    for(nodeID = 1; nodeID < nodesPerElement+1; nodeID++){
        nds1 = this->elementMap.GetElement(elementID,nodeID);// extracting node id
            // lets sweep the degrees of fredom in each node, ND
            for(id = 0; id < ND; id++){
                //ie = (iin-1)*ND + id;
                //ig = (this->elementMap.GetElement(elementID,iin))*ND + id;
                exyVec1.push_back(this->gXYvec.GetElement(nds1,id)); // id-th coordinate
                euVec1.push_back(this->gu.GetElement(ND*(nds1)+id,0)); // displancement on the id-th axis
                // delete vectors
            } // sweeping degrees of freedom in each node
            //
    } // sweeping nodes of each element
    //
    qbMatrix2<T> exy(NF,ND,exyVec1);
    //cout << endl;
    //qbMatrix2<T> eu2(NEV,1,euVec1);
    eu.SetToZero();
    for(int ii = 0; ii < NEV; ii++){
    eu.SetElement(ii,0,euVec1.at(ii));
    }    
    //cout << endl;
    //eu.PrintMatrix();
    //cout << endl;
    //gstrain.PrintMatrix();
    //////////////////////////// lets loop over the gauss points
    for(int ig = 0; ig < ngp; ig++){
        lgp++;
        // lets compute local Cartesian Strains
        // lets compute the B matrix without aking for the Jacobian determinant
        B.SetToZero();
        B = this-> Bmtx_quad(pgp.GetElement(ig,0), pgp.GetElement(ig,1), exyVec1);
        //qbMatrix2<T> B = this->Bmtx_quad(pgp.GetElement(ig,0), pgp.GetElement(ig,1), exyVec1,detJ);
        //B.PrintMatrix();
        estrain = B*eu; // computing strain at gauss point
        //cout << endl;
        //estrain.PrintMatrix();
        //cout << endl;
        // lets compute local Cartesian Stresses
        this->Eval_D_vec(pgp.GetElement(ig,0), pgp.GetElement(ig,1));
        estress = D*estrain; // computing stress at gauss point
        //computing principal strains at gauss point
        quad<T> quad1(pgp.GetElement(ig,0), pgp.GetElement(ig,1));
        quad1.Eval_principalStress(estress);
        pstress.SetElement(0,0,quad1.get_pr1());
        pstress.SetElement(1,0,quad1.get_pr2());
        pstress.SetElement(2,0,quad1.get_alpha());
        quad1.Eval_local_sf(pgp.GetElement(ig,0), pgp.GetElement(ig,1));
        qbMatrix2<T> shp(1,NF,quad1.get_sF_vec());
        pXYvec = shp*exy;
        //pXYvec.PrintMatrix();
        //cout << endl;
        //
        // lets populate global Cartesian Strains and Stresses
        //for(nodeID = 1; nodeID < nodesPerElement+1; nodeID++){
            nds1 = this->elementMap.GetElement(elementID,ig+1);// extracting node id
                // lets sweep the degrees of fredom in each node, NS number of stress components
                for(id = 0; id < NS; id++){
                    localii = (ig)*NS + id;
                    globalii = (nds1)*NS + id;
                    //cout << ig << " " << nds1 << " "<< localii << " " << globalii << " " << estrain.GetElement(id,0) << endl;
                    gstrain.SetElement(globalii,0,estrain.GetElement(id,0));
                    gstress.SetElement(globalii,0,estress.GetElement(id,0));
                    gpstress.SetElement(globalii,0,pstress.GetElement(id,0));
                    if(id<ND) gpXYvec.SetElement(nds1,id,pXYvec.GetElement(0,id)); // Saving global gauss points
                    // delete vectors
                } // sweeping degrees of freedom in each node
                //
            cout << elementID+1 << setw(9) << nds1+1 << setw(8) << lgp << setw(15) << pXYvec.GetElement(0,0) << setw(15) << pXYvec.GetElement(0,1) << setw(10) << estress.GetElement(0,0) << setw(10) << estress.GetElement(1,0) << setw(8) << estress.GetElement(2,0) << setw(12) << pstress.GetElement(0,0) << setw(16) << pstress.GetElement(1,0)<< endl;
        //} // sweeping nodes of each element
        //
    //
    //
    } // sweeping gauss points
    //
    //
    // Deleting local elemnt strains and stress
    while (euVec1.size() > 0) {
        exyVec1.pop_back();
        euVec1.pop_back();
    }
    //
} // sweeping elements
//
//

//
//
}


// Function for evaluating matrix D
template <class T> void quadModel<T>::Eval_D_vec(T E1, T nu1)
{
T coef;
coef = E1/(1-nu1*nu1);
//
D_vec.push_back(coef);          //11
D_vec.push_back(coef*nu1);       //12
D_vec.push_back(0);             //13
D_vec.push_back(coef*nu1);       //21
D_vec.push_back(coef);          //22
D_vec.push_back(0);             //23
D_vec.push_back(0);             //31
D_vec.push_back(0);             //32
D_vec.push_back(coef*(1-nu1)/2); //33
}



template <class T> qbMatrix2<T> quadModel<T>::Ematx_quad(int nev, int ng, vector<T> xye, qbMatrix2<T> D, T t)
{
    /*
    nev = nn*nd : number of variables in an element
    nn : number of nodes per element
    nd : number of coordinate dimensions / degrees of freedom
    ns : number of stress components
    ng : number of Gaussian integration points
    t : element thickness
    */
    //
    T xi, eta;
    T detJ; // Jacobian matrix determinant
    //
    vector<T> wgp_vec; // vector of integration weights
    T c_val = 0.577350269189626; // CHECK accuracy
    vector<T> pgp_vec; // vector of Gauss Coordinates
    //
    T dv; // transformation factor for integration
    //
    // initializizng stiffness matrix
    qbMatrix2<T> Ke(nev,nev);
    //
    //
    // Set gauss point positions and weights
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(-c_val);    
    //
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(-c_val);
    //
    pgp_vec.push_back(c_val);
    pgp_vec.push_back(c_val);
    //
    pgp_vec.push_back(-c_val);
    pgp_vec.push_back(c_val);
    qbMatrix2<T> pgp(NF,ND,pgp_vec);
    //cout << "Printing Gauss Integration Coordinates: ";
    //pgp.PrintMatrix(); 
    //
    //
    for(int ii=0;ii<ng;ii++){
        wgp_vec.push_back(1.0);
    }
    //qbVector<T> wgp(wgp_vec); // vector of integration weights
    //
    //
    Ke.SetToZero();
    //cout << "Printing Ke stiff matrix: ";
    //Ke.PrintMatrix(); 
    //
    ///// Gaussina Integration Loop
    for(int ig=0;ig<ng;ig++){
        // selecting xi and eta values
        //xi = pgp.GetElement(ig,0);
        //eta = pgp.GetElement(ig,1);
        //cout << "" <<xi <<" "<< eta<<" "<< wgp.GetElement(ig)<<endl;
        qbMatrix2<T> B = this-> Bmtx_quad(pgp.GetElement(ig,0), pgp.GetElement(ig,1), xye,detJ);
        qbMatrix2<T> B_t = B.Transpose();
        //B.PrintMatrix();
        //cout << detJ << endl;
        //cout << wgp_vec[ig] << endl;
        dv = wgp_vec[ig]*detJ*t;
        //T dv = wgp.GetElement(ig)*detJ*t;
        //Ke = Ke + B_t*D*B*dv;
        Ke = Ke + B.Transpose()*D*B*dv;
        /*
        % Gaussian integration loop
        for ig=1:ng
            xi=pgp(ig,1) ; eta=pgp(ig,2);
            [shp,sdv]=Shapf(xi,eta); %#ok<ASGLU>
            [B,detJ]=JBMat(xyel,sdv);
            dv=wgp(ig)*detJ*t;
            Ke=Ke+B'*D*B*dv;
        end
        */
       //
    }
    //
    //
    return Ke;
}

//returns the B matrix for a quad linear element, returning Jacobian determinant
template <class T> qbMatrix2<T> quadModel<T>::Bmtx_quad(T xi, T eta, vector<T> xye, T &J1Det)
{
    //T J1Det;
    qbMatrix2<T> J1,invJ1;
    qbMatrix2<T> sd1;
    qbMatrix2<T> preB(NS,NEV);
    qbMatrix2<T> Bm2D(NS,NEV);
    //
    //
    // computing shape functions, 1x4
    quad<T> quad1(xi,eta);
    //
    // computing local shape functions derivatives, 2x4
    //cout << "Printing Shape Function Derivatives Matrix dsF1:";
    quad1.Eval_local_sd(xi,eta);
    qbMatrix2<T> dsF1(ND,NF,quad1.get_dsF_vec());
    //dsF1.PrintMatrix();
    //
    // computing Jacobian, 2x2
    //cout << "Printing global Jacobian Matrix J1:";
    qbMatrix2<T> XY1(NF,ND,xye);
    J1 = dsF1*XY1;
    //J1.PrintMatrix();
    //
    // computing Jacobian determinat
    J1Det = J1.Determinant();
    //cout << "Printing global Jacobian Matrix Determinant: "<< J1Det <<"";
    invJ1 = J1; // initializing inverse matrix
    if(invJ1.Inverse()){
        //cout << "Printing global Invrser Jacobian Matrix invJ1: "; // Inverse Jacobian, 2x2 matrix
        //invJ1.PrintMatrix();
        //cout << "Printing global Cartesian Derivatives sd1:";     // Global Cartesain derivatives, 2x4 matrix
        sd1 = invJ1*dsF1;
        //sd1.PrintMatrix();
        //cout << "Creating B matrix with custom class  qbMatrix2 :"; // B matrix, 3x8
        Bm2D.SetToZero();
        //cout << endl;
        //Bm2D.PrintMatrix();
        //sd1.PrintMatrix();
        //
        Bm2D.SetElement(0,0,sd1.GetElement(0,0));
        Bm2D.SetElement(0,1,0);
        Bm2D.SetElement(0,2,sd1.GetElement(0,1));
        Bm2D.SetElement(0,3,0);
        Bm2D.SetElement(0,4,sd1.GetElement(0,2));
        Bm2D.SetElement(0,5,0);
        Bm2D.SetElement(0,6,sd1.GetElement(0,3));
        Bm2D.SetElement(0,7,0);
        //
        Bm2D.SetElement(1,0,0);
        Bm2D.SetElement(1,1,sd1.GetElement(1,0));
        Bm2D.SetElement(1,2,0);
        Bm2D.SetElement(1,3,sd1.GetElement(1,1));
        Bm2D.SetElement(1,4,0);
        Bm2D.SetElement(1,5,sd1.GetElement(1,2));
        Bm2D.SetElement(1,6,0);
        Bm2D.SetElement(1,7,sd1.GetElement(1,3));
        //
        Bm2D.SetElement(2,0,sd1.GetElement(1,0));
        Bm2D.SetElement(2,1,sd1.GetElement(0,0));
        Bm2D.SetElement(2,2,sd1.GetElement(1,1));
        Bm2D.SetElement(2,3,sd1.GetElement(0,1));
        Bm2D.SetElement(2,4,sd1.GetElement(1,2));
        Bm2D.SetElement(2,5,sd1.GetElement(0,2));
        Bm2D.SetElement(2,6,sd1.GetElement(1,3));
        Bm2D.SetElement(2,7,sd1.GetElement(0,3));
        //
        //Bm2D.PrintMatrix();  
        //
        /*
        quad1.Eval_Jm2D(eXYvec1);
        //quad1.print_Jm2D();
        //qbMatrix2<T> J2(ND,ND,quad1.get_Jm2D_vec());
        //J2.PrintMatrix();
        // 
        quad1.Eval_inv_Jm2D();
        //quad1.print_invJm2D();
        quad1.Eval_global_sd();
        //quad1.print_global_sd();
        quad1.Eval_Bm2D();
        //quad1.print_bm2D();
        qbMatrix2<T> Bm2D_b(NS,NEV,quad1.get_B_vec());
        //Bm2D_b.PrintMatrix();  
        */
        //
        //
        return Bm2D;
    }
    else{
        throw std::invalid_argument("Cannot compute the inverse of a matrix .");
    }    
}
///////////////////////////////////////// Function overloading, no Jacobian Determinant
template <class T> qbMatrix2<T> quadModel<T>::Bmtx_quad(T xi, T eta, vector<T> xye)
{
    //T J1Det;
    qbMatrix2<T> J1,invJ1;
    qbMatrix2<T> sd1;
    qbMatrix2<T> preB(NS,NEV);
    qbMatrix2<T> Bm2D(NS,NEV);
    Bm2D.SetToZero();
    //
    //
    // computing shape functions, 1x4
    quad<T> quad1(xi,eta);
    //

    // computing local shape functions derivatives, 2x4
    //cout << "Printing Shape Function Derivatives Matrix dsF1:";
    quad1.Eval_local_sd(xi,eta);
    qbMatrix2<T> dsF1(ND,NF,quad1.get_dsF_vec());
    
    //dsF1.PrintMatrix();
    //
    // computing Jacobian, 2x2
    //cout << "Printing global Jacobian Matrix J1:";
    //for(auto ii = xye.begin(); ii != xye.end(); ii++){
    //    cout << *ii << " "<< distance(xye.begin(),ii) << endl;
    //}
    qbMatrix2<T> XY1(NF,ND);
    //XY1.SetToZero();
    //for(auto ii = xye.begin(); ii != xye.end(); ii++){
    //    int jj = distance(xye.begin(),ii);
    //    XY1.SetElement(jj%NF-jj%ND,jj%ND,*ii);
    //    cout << (jj-jj%NF)<< " "<<jj%ND << " " << XY1.GetElement(jj%NF,jj%ND) << " " << *ii << " "<< distance(xye.begin(),ii) << endl;
    //}
    auto kk = xye.begin();
    for(int ii=0;ii<NF;ii++){
        for(int jj=0;jj<ND;jj++){
            if(kk<xye.end()){
            //XY1.SetElement(ii,jj,xye.at((jj*NF)+ii));
            XY1.SetElement(ii,jj,*kk);
            kk++;
            }
        }
    }
    //XY1.PrintMatrix();
    //cout<<"\n";
    J1 = dsF1*XY1;
    //J1.PrintMatrix();
    //

    //cout << "Printing global Jacobian Matrix Determinant: "<< J1Det <<"";
    invJ1 = J1; // initializing inverse matrix
    //invJ1.PrintMatrix();
    
    if(invJ1.Inverse()){
        //cout << "Printing global Invrser Jacobian Matrix invJ1: "; // Inverse Jacobian, 2x2 matrix
        //invJ1.PrintMatrix();
        //cout << "Printing global Cartesian Derivatives sd1:";     // Global Cartesain derivatives, 2x4 matrix
        sd1 = invJ1*dsF1;
        //sd1.PrintMatrix();
        //cout << "Creating B matrix with custom class  qbMatrix2 :"; // B matrix, 3x8
        Bm2D.SetToZero();
        //cout << endl;
        //Bm2D.PrintMatrix();
        //sd1.PrintMatrix();
        //
        Bm2D.SetElement(0,0,sd1.GetElement(0,0));
        Bm2D.SetElement(0,1,0);
        Bm2D.SetElement(0,2,sd1.GetElement(0,1));
        Bm2D.SetElement(0,3,0);
        Bm2D.SetElement(0,4,sd1.GetElement(0,2));
        Bm2D.SetElement(0,5,0);
        Bm2D.SetElement(0,6,sd1.GetElement(0,3));
        Bm2D.SetElement(0,7,0);
        //
        Bm2D.SetElement(1,0,0);
        Bm2D.SetElement(1,1,sd1.GetElement(1,0));
        Bm2D.SetElement(1,2,0);
        Bm2D.SetElement(1,3,sd1.GetElement(1,1));
        Bm2D.SetElement(1,4,0);
        Bm2D.SetElement(1,5,sd1.GetElement(1,2));
        Bm2D.SetElement(1,6,0);
        Bm2D.SetElement(1,7,sd1.GetElement(1,3));
        //
        Bm2D.SetElement(2,0,sd1.GetElement(1,0));
        Bm2D.SetElement(2,1,sd1.GetElement(0,0));
        Bm2D.SetElement(2,2,sd1.GetElement(1,1));
        Bm2D.SetElement(2,3,sd1.GetElement(0,1));
        Bm2D.SetElement(2,4,sd1.GetElement(1,2));
        Bm2D.SetElement(2,5,sd1.GetElement(0,2));
        Bm2D.SetElement(2,6,sd1.GetElement(1,3));
        Bm2D.SetElement(2,7,sd1.GetElement(0,3));
        //
        //Bm2D.PrintMatrix();  
        //
        //
        //
        return Bm2D;
    }
    else{
        throw std::invalid_argument("Cannot compute the inverse of a matrix .");
    }    
}

// Define error codes.
constexpr int QBLINSOLVE_NOUNIQUESOLUTION = -1;
constexpr int QBLINSOLVE_NOSOLUTIONS = -2;

// The qbLinSolve function.
template <class T> int quadModel<T>::qbLinSolve(const qbMatrix2<T> &aMatrix, const qbVector<T> &bVector, qbVector<T> &resultVec)
{
	// Make a copy of the input matrix, aMatrix.
	// We will use this to create the augmented matrix, so we have
	// to make a copy.
	qbMatrix2<T> inputMatrix = aMatrix;
	
	// Compute the rank of the original matrix.
	int originalRank = inputMatrix.Rank();

	/* Combine inputMatrix and bVector together into a single matrix,
		ready for using Gaussian elimination to reduce to 
		row-echelon form. */
	
	// Extract data from bVector.
	int numDims = bVector.GetNumDims();
	std::vector<T> bVecData;
	for (int i=0; i<numDims; ++i)
		bVecData.push_back(bVector.GetElement(i));
		
	// Use this to create a qbMatrix2 object with the same data (nx1).
	qbMatrix2<T> bMatrix(numDims, 1, bVecData);
	
	// Combine the two matrices together.
	inputMatrix.Join(bMatrix);
	
	/* Use Gaussian elmination to convert to row-echelon form. */
	qbMatrix2<T> rowEchelonMatrix = inputMatrix.RowEchelon();
	
	/* Comute the rank of the augmented matrix.
		Note that we do this after performing Gaussian elimination to
		reduce the matrix to row echelon form so that if this was 
		successful, there is no need to repeat this operation twice. */
	int augmentedRank = rowEchelonMatrix.Rank();
	
	/* ********************************************************************* 
		Test the two ranks to determine the nature of the system we
		are dealing with. The conditions are as follows:
		
		n = number of rows.
		
		1) originalRank = augmentedRank = n	=> A unique solution exists.
		2) originalRank = augmentedRank < n	=> An infinite number of solutions exist.
		3) originalRank < augmentedRank			=> No solutions exist.  
		********************************************************************* */
	if ((originalRank == augmentedRank) && (originalRank < inputMatrix.GetNumRows()))
	{
		return QBLINSOLVE_NOUNIQUESOLUTION;
	}
	else if (originalRank < augmentedRank)
	{
		return QBLINSOLVE_NOSOLUTIONS;
	}
	else
	{
		/* Create a qbVector object to store the output. Initially we will
			populate this with the data from bVecData, but we are going to modify
			the elements as we compute them. */
		qbVector<T> output(bVecData);
		
		// Now use back-substitution to compute the result.
		int numRows = rowEchelonMatrix.GetNumRows();
		int numCols = rowEchelonMatrix.GetNumCols();
		int startRow = numRows-1;
		
		// Loop over the rows, in reverse order.
		for (int i=startRow; i>=0; --i)
		{
			// Extract the currentResult for this row.
			T currentResult = rowEchelonMatrix.GetElement(i, numCols-1);
	
			// Compute the cumulative sum.
			T cumulativeSum = static_cast<T>(0.0);
			for (int j=i+1; j<numRows; ++j)
			{
				cumulativeSum += (rowEchelonMatrix.GetElement(i,j) * output.GetElement(j));
			}
			
			// Compute the answer.
			T finalAnswer = (currentResult - cumulativeSum) / rowEchelonMatrix.GetElement(i,i);
			
			// And store.
			output.SetElement(i, finalAnswer);
			
		}
		
		// Return the output.
		resultVec = output;
		
	}
		
	return 1;	
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//deconstructor
template <class T> quadModel<T>::quadModel::~quadModel()
{

}







template class quad<long double>;
template class quad<double>;
template class quad<int>;

template class quadModel<long double>;
template class quadModel<double>;
template class quadModel<int>;
