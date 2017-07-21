

#include <vector>
#include <iostream>
#include <math.h>


#include "vecmanip.h"
#include "py_solution.h"
#include "potential.h" 
using namespace std;

int main()
{

	double rho = 0.1;
	double rm = 10.;
	int N = 5000;
	int Nit = 100;
	double a = .9;
	

	double eps = 1.;
	double beta = 1.;
	double sigma = 1.;
	double rcu_min = 0.1;
	double rcu_max = pow(2.,1./6)*sigma;

	double rm2 = 2*rm;
	int N2 = 2*N;
	double dr = rm/(N-1);

	vector<double> e_minbu = getvec_exp_minbu(rm2,N2,eps,beta,sigma,rcu_min,rcu_max);

	vector<double> h0 = get_h0(N2,dr);
	
	vector<double> g0 = get_h0(N2,dr);
	h2g(g0,e_minbu,dr,rcu_min);
	write_vec(g0,"g0.dat");

	get_h(h0,e_minbu,rho,dr,N,Nit,a);

	h2g(h0,e_minbu,dr,rcu_min);

	vector<double> r(N2);
	for(int ri=0;ri<N2;++ri)
		r[ri] = ri*dr;


	write_vec(r,"r.dat");
	write_vec(h0,"gr.dat");
	return 0;
}


