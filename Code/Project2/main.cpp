#include <iostream>
#include <fstream>
#include <armadillo>
#include "functions.h"
#include "test.h"

using namespace std;
using namespace arma;


int main()
{


    cout<<"running test"<<endl;
    testing();
    cout<<"Finnished test"<<endl;
    int n = 3;
    int a = 3;
    int d = 2;
    ofstream file, outfile;
    outfile.open("/Users/andreas/Documents/Fys4150/Project-2/Code/Project2/mytext.txt", ios_base::trunc);

    mat test = tridiagonal(n,a,d);
    outfile << test<< endl;
    vec eigenvalues = analytic_eigenvalues(n,a,d);
    outfile << eigenvalues;
    mat R = mat(n,n,fill::eye);
    //R = jacobi_method(test,R,n);
    outfile << R;

    outfile.close();
    cout<<"Finnished"<<endl;
    /**/
    return 0;
}
