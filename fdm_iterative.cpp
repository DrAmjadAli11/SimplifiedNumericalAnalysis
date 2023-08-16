#include<iostream>
#include<cmath>
#include<iomanip>
using namespace std;

#define a		1.0
#define b		2.0
#define alpha	1.0
#define beta	2.0
#define m		9		//number of interior nodes
#define N		200		//maximum number of iterations
#define TOL		0.0000001

#define p(x) -2.0/(x)
#define q(x) 2.0/(x*x)
#define r(x) sin(log(x))/(x*x) 

void egs(double [], double [] , double [] , 
			double [], double [], double h) ;


int main()
{  
	double h , err , sum ;
	double x[m+2] , z[m+2] , zp[m+2] ;
	double B[m+1] , D[m+1] , U[m+1], L[m+1] ;
	int i , k ;
	
	h = (b-a)/(m+1) ;
	x[0] = a ;
	x[m+1] = b ;

	for(i=1 ; i<=m ; i++)
	{
		x[i] = x[i-1]+h ;

		cout<< "\tnodes" << x[i] << endl ;
	}
	
	z[0] = alpha ;
	z[m+1] = beta ;
	
	for (i=1 ; i<=m ;i++)
		B[i] = (-h)*(h)*(r(x[i]));
	
	for (i=1 ; i<=m ; i++)
		D[i] = 2+((h)*(h)*(q(x[i])));

	for (i=1 ; i<=m ; i++)
		U[i] = -1.0 + ((h*0.5)*(p(x[i])));
   
	for (i=1 ; i<=m ; i++)
		L[i] = -1.0-((h*0.5)*(p(x[i]))) ;
    
    for(i=1 ; i<=m ; i++)
	   	z[i] = 0.0 ;

	k = 0 ;
	cout<< setw(4) << k << ": " ;
	cout<< "z= " ;
	cout<< setprecision(2) << setw(3) << z[0] << " " ;
	for (i=1 ; i<=m ; i++)
	{
		cout<< setprecision(8) << setw(9) << z[i] << " " ;
	}
	cout<< setprecision(2) << setw(3) << z[m+1] << " " ;
	cout << endl ;
	
	egs(z, B, D, U, L, h) ;
	
	return 0 ;
}

//  ---- User-defined function for an efficiant Gauss-Seidel method ----  //

void egs(double z[m+2], double B[m+1] , double D[m+1] , 
			double U[m+1], double L[m+1], double h)
{
	int i, k ;
	double sum , err ;
	double zp[m+2] ;

	for (k=1 ; k<=N ; k++)        // Iterations loop
	{
   		// Making a copy of the solution vector before updating it
		for (i=1 ; i<=m ; i++)   zp[i] = z[i] ;
   
   		// Updating the solution vector
   		for (i=1 ; i<=m ; i++)
		{
			z[i] = (B[i]-(L[i]*z[i-1])-(U[i]*z[i+1]))/D[i] ;
		}
		
		// Printing the latest solution vector
		cout<< setw(4) << k << ": " ;
		cout<< "z= " ;
		cout<< setprecision(2) << setw(3) << z[0] << " " ;
		for (i=1 ; i<=m ; i++)
		{
			cout<< setprecision(8) << setw(9) << z[i] << " " ;
		}
		cout<< setprecision(2) << setw(3) << z[m+1] << " " ;
		cout << endl ;
		
		// Finding the error as the L2-norm
		sum = 0.0 ;
		for (i=1 ; i<=m ; i++)
			sum = sum + (z[i]-zp[i])*(z[i]-zp[i])/(z[i]*z[i]) ;
		err = sqrt(sum) ;
		
		// Testing the convergence
		if (err<TOL) break ;
		
	}
}
