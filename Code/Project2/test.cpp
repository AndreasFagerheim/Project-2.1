#include <iostream>
#include <fstream>
#include <armadillo>
#include "functions.h"
#include "test.h"
#include <cmath>


using namespace std;
using namespace arma;

int testing(){
    cout <<"test initiated!"<<endl;

    // test if maxoffdiag finds the largets value
    int n = 5;
    double eps = pow(10,-8);
    double maxValue = -51.9;
    mat temp = mat(n,n,fill::eye);

    temp(2,3) = 5;
    temp(3,2) = 5;
    temp(1,4) = 51.2;
    temp(4,1) = 51.2;
    temp(3,0) = -51.9;
    temp(0,3) = -51.9;
    int k = 0;
    int l = 0;
    double maxfound = maxoffdiag(temp,&k,&l,n);
    cout<<"Testing max value"<<endl;
    if (maxfound == fabs(maxValue)){
        cout<<"Success! Max value = expected max value"<<endl;
    }else{
        cout<<"ERROR!!! (Not the right max value)"<<endl;
    }


    // testing eignevalues of implemeted algorithm vs analytic vs armadillo
    int a = 2;
    int d = 3;
    mat testMatrix = tridiagonal(n, a,d);
    // eigenvalues using armadillo
    vec eig_armadillo = sort(eig_sym(testMatrix));
    // eigenvalues grom analytic solution
    vec eig_analytic = sort(analytic_eigenvalues(n,a,d));
    // eigenvalues from jacobis method
    vec eig_jacobi = 0;

    // analytic vs armadillo method
    for(int i = 0;i<n; i++){
        if(fabs(eig_analytic[i]-eig_armadillo[i]) < eps){
            cout<<"Sucsess!! Analytic eigenvalues equal to armadillos";
        }else{
            cout<<"ERROR!! Differnet eigenvalues for analytic and armadillos method";
        }
    }

    return 0;
}
