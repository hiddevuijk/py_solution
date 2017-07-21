#ifndef GUARD_potential_h
#define GUARD_potential_h

#include <math.h>
#include <vector>



double u(double r,double eps, double sigma)
{
	double sr6 = pow(sigma/r,6.);
	return 4*eps*(sr6*sr6-sr6);
}

// if rcu_max <= 0, domain u(r) is (rcu_min,inf)
double exp_minbu(double r,double eps,double beta, double sigma, 
	double rcu_min, double rcu_max)
{
	double u_rcu = rcu_min<rcu_max ? u(rcu_max,eps,sigma) : 0.;
	if(r<rcu_min) return 0.;
	if(r>rcu_max and rcu_max>rcu_min) return 1.;
	return exp(-1*beta*u(r,eps,sigma) +  beta*u_rcu);
}

std::vector<double> getvec_exp_minbu(double rmax, int N,
	double eps,double beta, double sigma, double rcu_min, double rcu_max)
{
	std::vector<double> temp(N);
	double dr = rmax/(N-1);
	for(int i=0;i<N;++i)
		temp[i] = exp_minbu(i*dr,eps,beta,sigma,rcu_min,rcu_max);
	return temp;
}
 




#endif
