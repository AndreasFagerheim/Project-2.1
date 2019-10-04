#include "functions.h"

#include <iostream>
#include <cmath>
#include <armadillo>


using namespace std;
using namespace arma;

/*
 *  A = input matrix (n x n)
 *  R = empty matrix, eigenvectors (n x n)
 *  n = dimentino of matrix
 *
*/

//Function creating a tridiagonal matrix with diagonal (d) and second diagonal (a)

mat tridiagonal(int n, double a, double d){
    mat A = mat(n,n,fill::zeros);
    for(int i=0;i<n;i++){
        for(int j = 0;j<n;j++){
            if(i==j){
                A(i,j) = d;
            }else if(i==1+j){
                A(i,j) = a;
            }else if(i == j-1){
                A(i,j) = a;
            }
        }
    }
    return A;
}
// returns an array with eiegnevalues given from analytical solutions
vec analytic_eigenvalues(int n, double a, double d){
    vec v = vec(n);
    double pi = 4.0*atan(1);
    for(int i = 0; i<n;i++){
        v(i) = d+2 * a * cos((i * pi)/(n+1)); // given in the project sheet
    }
    return v;
}

// Find biggest matrix element in the upper half
double maxoffdiag(mat &A, int * k, int * l, int n){
    double max = 0.0;
    for (int i = 0; i<n;i++){
        for(int j= i+1;j<n;j++){
            if(fabs(A(i,j))> max){
                max = fabs(A(i,j));
               *k = i;
               *l = j;
            }
        }
    }
    return max;
}


//Find values of sin and cos

void rotate(mat &A, mat &R, int k, int l, int n){

    double s,c;
    if(A(k,l) != 0.0){
        double t, tau;
        tau = (A(l,l)-A(k,k))/(2.0*A(k,l));
        if(tau > 0){
            t = 1.0/(tau + sqrt(1.0 + tau*tau));
        } else{
            t = -1.0/(-tau +sqrt(1.0 + tau*tau));
        }
        c = 1.0/sqrt(1+t * t);
        s = c * t;
    }else{
        c = 1.0;
        s = 0.0;
    }
    double a_kk,a_ll,a_ik,a_il,r_ik,r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);

    // changing the matrix elements with indices k and l
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + s*s*a_ll;
    A(k,l) = 0.0; // set them manually to zero
    A(l,k) = 0.0; // set them manually to zero

    // then change remaining elements
    for(int i = 0;i<n;i++){
        if(i!=k && i!=l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        // Eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
   // cout << A<<endl;
   // cout <<"R"<< R<<endl;
    return;

}
/*
// Find biggest matrix element in the upper half
void max2(mat &A, int &k, int & l, double epsilon, bool &achieved){
    double max = 0.0;
    int n = A.n_rows;
    for (int i = 0; i<n;i++){
        for(int j= 0;j<n;j++){
            if(abs(A(i,j)>abs(max))){
                max = abs(A(i,j));
                k = i;
                l = j;
            }
        }
    }
}
/**/
/*
// alternative way of rotate testing????
void rotation2(mat &A, double c, double s, int k , int l){
    double a_ik,a_il,a_ki,a_li;
    for(int i = 0; i<A.n_rows;i++){
        a_ki = A(k,i);
        a_li = A(l,i);
        A(k,i) = a_ki*c - a_li*s;
        A(l,i) = a_ki*s + a_li*c;
    }
    for(int i = 0; i<A.n_rows;i++){
        a_ik = A(i,k);
        a_il = A(i,l);
        A(i,k) = a_ik*c - a_il*s;
        A(i,l) = a_ik*s + a_il*c;
    }
}
vec jacobi_method2(mat A, double epsilon, int &counter){
    bool achieved = false;
    int k = 0;
    int l = 0;
    double t,tau,c,s;

}
/**/
mat jacobi_method (mat A, mat R, int n){
    // setting up eigenvector matrix
    for (int i = 0; i < n; i++){
        for (int j = 0; j<n; j++){
            if (i == j){
                R(i,j) = 1.0;
            }else{
                R(i,j) = 0.0;
            }
        }

    }
    cout <<R<<endl;
    cout <<A<<endl;

    int k,l;
    double epsilon = 1.0e-8;
    double max_number_iterations = (double) n * (double)  n *(double)  n;
    int iterations = 0;
    double max_offdiag = maxoffdiag(A,&k,&l,n);

    while (abs(max_offdiag)>epsilon && double (iterations) <max_number_iterations){
        max_offdiag = maxoffdiag(A,&k,&l,n);
        rotate(A, R, k, l, n);
        iterations++;
    }
    cout<<"Number of iterations:" <<iterations<<"\n";
    cout <<R<<endl;
    cout <<A<<endl;
}



