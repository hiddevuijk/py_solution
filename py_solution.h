#ifndef GUARD_py_solution_h
#define GUARD_py_solution_h

#include <vector>
#include <iostream>
#include <assert.h>

namespace PY {
	double pi2 = 2*acos(-1.);
}

double integrate_midpoint(const std::vector<double>& v,
	double dr, int il, int ir)
{
	if(il == ir) return 0.;
	double temp =0.; 
	for(int i=il;i<ir;++i)
		temp += v[i];
	return dr*temp;
}

void get_hnew(
	std::vector<double>& hnew,
	const std::vector<double>& hold,
	const std::vector<double>& e_minbu,
	std::vector<double>& hp,
	std::vector<double>& integrant,
	int N,double ds, double rho)
{
	int s_smr;
	double intsi;
	for(int ri=0;ri<N;++ri){
		// update integrant
		for(int si=0;si<N;++si) {
			s_smr = si < ri ? -1 : 1;
			intsi = hold[si+ri]*e_minbu[si+ri];
			intsi += s_smr*hold[s_smr*(si-ri)]*e_minbu[s_smr*(si-ri)];
			intsi -= 2*si*ds;
			integrant[si] =intsi*(1-e_minbu[si])*hold[si];
		}
		hp[ri] = 1-PY::pi2*rho*integrate_midpoint(integrant,ds,0,N);
		
	}
	hnew[0] = 0;
	for(int ri=1;ri<N;++ri) {
		//hnew[ri] = integrate_midpoint(hp,ds,0,ri);
		hnew[ri] = hnew[ri-1]+integrate_midpoint(hp,ds,ri-1,ri);
	}
}


void get_h(std::vector<double>& h0,
	const std::vector<double>& e_minbu,
	double rho, double dr, int N,
	int Nit, double a)
{
	std::vector<double> hnew(2*N);
	std::vector<double> hp(2*N);
	std::vector<double> integrant(2*N);
	for(int it=0;it<Nit;++it) {
		get_hnew(hnew,h0,e_minbu,hp,integrant,N,dr,rho);
		for(int ri=0;ri<N;++ri)
			h0[ri] = a*h0[ri] + (1-a)*hnew[ri];
	}

}

void h2g(std::vector<double>& h,
	const std::vector<double>& e_minbu,
	double dr,double rcu_min)
{
	for(int ri=0;ri<h.size();++ri) {
		if(ri*dr<rcu_min) h[ri] = 0.;
		else h[ri] *= e_minbu[ri]/(ri*dr);
	}
}	

std::vector<double> get_h0(int N2, double dr)
{
	std::vector<double> h0(N2);
	for(int ri=0;ri<N2;++ri)
		h0[ri] = ri*dr;
	return h0;
}












#endif
