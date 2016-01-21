/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_2d_e.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Mls_2d_e_exp : public Simulator < real, dim, Mls_2d_e_exp<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Matrix<real, 5, 1> vec5;
		typedef Eigen::Matrix<real, 5, 5> mat55;
		typedef Eigen::Triplet<real> tpl;
	public:
		Mls_2d_e_exp() {}
		~Mls_2d_e_exp() {}

		void step() {
			calInvMat();
			convect1();
			solPres();
			convect2();
			calInvMat();
			calPnd();
			makeFs();
			check();
			shift();
			//calInvMat(); //for sensor
		}

		void convect1() {
			/*standard*/
			/*
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] += para.dt * part->vel2[p];
			}
			*/
			//-------------------------------------------------------------------------//
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
		}

		void convect2() {
			/*standard*/
			/*
			std::vector<vec> dash(part->np);
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			dash[p] = -para.dt / para.rho * part->grad_mls_poly2d_d(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel1[p] = part->vel2[p] = part->vel2[p] + dash[p];
			part->pos[p] += para.dt * dash[p];
			}
			*/
			//------------------------------------------------------------------------//
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -(para.dt / para.rho) * part->grad(part->pres, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
				part->vel1[p] = part->vel2[p];
			}
			//------------------------------------------------------------------------//
			/*fourth-order Ronge Kuta*/
			/*
			std::vector<vec> x0(part->np);
			std::vector<vec> k2(part->np);
			std::vector<vec> k3(part->np);
			std::vector<vec> k4(part->np);
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			x0[p] = part->pos[p];
			part->pos[p] += 0.5* para.dt* part->vel1[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k2[p] = part->vel2[p] + 0.5* para.dt* (-1. / para.rho)* part->grad_mls_poly2d_d(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + 0.5* para.dt* k2[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k3[p] = part->vel2[p] + 0.5* para.dt* (-1. / para.rho)* part->grad_mls_poly2d_d(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + para.dt* k3[p];
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			k4[p] = part->vel2[p] + para.dt* (-1. / para.rho)* part->grad_mls_poly2d_d(part->pres, p);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->pos[p] = x0[p] + (para.dt / 6.)* (part->vel1[p] + 2.*k2[p] + 2.*k3[p] + k4[p]);
			}
			for (unsigned p = 0; p < part->np; p++) {
			if (part->type[p] != FLUID) continue;
			part->vel2[p] += -(para.dt / para.rho) * part->grad_mls_poly2d_d(part->pres, p);
			part->vel1[p] = part->vel2[p];
			}
			*/
			//------------------------------------------------------------------------//
		}

		void solPres() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				const real gamma = 3.;
				const real cs = 40.;
				const real coef = para.rho *cs*cs / gamma;
				const real ratio = 1. - para.dt* part->div(part->vel2, p);
				part->pres[p] = coef * (pow(ratio, gamma) - 1.);
			}
		}

		void init_() {
			part = new Particle_2d_e<real, dim>(para.k, para.beta);
			part->clean();
			*part << "Geo.in";
			part->buildCell();
			part->b2b();
			part->b2norm();
			//part->updateTeam();
			part->initInvMat();
		}

	public:
		Particle_2d_e<real, dim>*  part;

	private:
		std::vector<tpl> coef;

	};

}