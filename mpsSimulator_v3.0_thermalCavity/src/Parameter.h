/*
*/
#pragma once
#include "Header.h"

namespace SIM {

	template <typename R, unsigned D>
	class Parameter {
		typedef Eigen::Matrix<R,D,1> Vec;
	public:
		Parameter() { init(); }
		~Parameter() {}

	public:
		R k;		//radius k*dp
		R rho;	//density
		R niu;	//kinetic viscosity
		R dtMax;	//max time step
		R cfl;	//cfl num
		R tt;	//total time
		R eps;	//eps
		R beta;	//beta: free surface
		R alpha;	//particle shifting
		R c;		//particle shifting
		
		Vec g;		//gravity
		R dt;	//current time step
		R umax;	//current umax
	private:
		void init() {
			dt = 1.;
			g = Vec::Zero();
		}
	};

}