/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include "Header.h"
#include "Parameter.h"
#include "Particle.h"
#include "MatSolver.h"
#include "Shifter.h"
#include "Sensor.h"

namespace SIM {

	template <typename R, unsigned D, typename Derived>
	class Simulator {
		typedef Eigen::Matrix<R,D,1> Vec;
		typedef Eigen::Matrix<R,D,D> Mat;
		typedef Eigen::Triplet<R> Tpl;
	public:
		Simulator() { timeStep = 0; }
		~Simulator() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void operator >> (const std::string& str) const {
			saveData(str);
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " No file Para. found ! " << std::endl;
			file >> para.k >> para.Pr >> para.Ra >> para.cfl >> para.dtMax >> para.tt >> para.eps >> para.alpha >> para.beta;
			std::cout << " Effective radius (times of dp)   : " << para.k << std::endl;
			std::cout << " Prandtl number                   : " << para.Pr << std::endl;
			std::cout << " Rayleigh number                  : " << para.Ra << std::endl;
			std::cout << " CFL number                       : " << para.cfl << std::endl;
			std::cout << " Maximum time step (1)            : " << para.dtMax << std::endl;
			std::cout << " Total time (1)                   : " << para.tt << std::endl;
			std::cout << " EPS                              : " << para.eps << std::endl;
			std::cout << " Arbitrary parameter Alpha        : " << para.alpha << std::endl;
			std::cout << " Arbitrary parameter Beta         : " << para.beta << std::endl;
			std::cout << " Reading Para. done " << std::endl;
			file.close();
		}

		void init() {
			derived().init_();
			mSol = new MatSolver<R, D>(unsigned(derived().part->np), para.eps);
			std::cout << " Particle number : " << derived().part->np << std::endl;
			*(derived().sen) << "Sensor.in";
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			timeStep = int(derived().part->ct / para.dt);
		}

		void mainLoop() {
			auto* const part = derived().part;
			while (part->ct <= para.tt) {
				std::cout << " step ----------------------------------> " << timeStep << std::endl;
				R tmp = cfl();
				para.dt = tmp < para.dtMax ? tmp : para.dtMax;
				part->updateCell();
				derived().step();
				part->ct += para.dt;	timeStep++;
				std::cout << " time --------> " << part->ct << std::endl;
				std::cout << " dt ----------> " << para.dt << std::endl;
			}
			saveData();
			sensorOut();
		}

		R stepGL() {
			static const R Re = R(1.) / R(para.Pr);
			static int counter = 0;
			static int maxLoop = 1;
			static R minDt = 0.1* cfl();
			auto* const part = derived().part;
			if (part->ct > para.tt) {
				saveData();
				sensorOut();
				std::exit(0);
			}
			std::cout << " step ----------------------------------> " << timeStep << std::endl;
			R tmp = cfl();
			para.dt = tmp < para.dtMax ? tmp : para.dtMax;
			part->updateCell();
			derived().step();
			part->ct += para.dt;	timeStep++;
			std::cout << " time --------> " << part->ct << std::endl;
			std::cout << " dt ----------> " << para.dt << std::endl;
			counter++;
			maxLoop = static_cast<int>(1.* Re / para.dt);
			if (counter > maxLoop && para.dt > minDt) {
				counter = 0;
				para.cfl /= R(2);
			}
			return part->ct;
		}

		void sensorOut() {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().sen) >> convert.str();
		}

		void profileOut() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			derived().sen->profile(rt, pf);
		}

		void profileOut_avgVel2() {
			auto* const part = derived().part;
			static std::string pf = "profile";
			auto sum = 0.;
			auto count = 0;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					sum += part->vel2[p].squaredNorm();
					count++;
				}
			}
			sum = sum / count;
			std::ofstream file("./out/" + pf + ".out", std::ofstream::app);
			file << std::setprecision(6) << std::scientific << timeStep*para.cfl << " "
				<< std::setprecision(6) << std::scientific << sum
				<< std::endl;
			file.close();
			std::cout << " Writing profile. done " << std::endl;
		}

		void saveData() const {
			static int i = 0;
			std::ostringstream convert;
			convert << i++;
			*(derived().part) >> ("./out/" + convert.str() + ".out");
		}
		void saveData(const std::string& str) const {
			*(derived().part) >> ("./out/" + str + ".out");
		}

		void fina() {}

	public:
		Parameter<R,D> para;
		MatSolver<R,D>* mSol;
		Shifter<R,D> shi;

	protected:
		void step() {}
		void convect() {}
		void visTerm_e() {}
		void visTerm_i() {}
		void presTerm_e() {}
		void presTerm_i() {}
#if LEGACY
		void makeDirichlet_p_op() {
			auto* const part = derived().part;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				part->fs[p] = 1;
				break;
			}
		}

		void makeDirichlet_p_avg() {
			auto* const part = derived().part;
			Eigen::SparseMatrix<R> d(part->np, part->np);
			std::vector<Tpl> coef;
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				for (int q = 0; q<int(part->np); q++) {
					if (part->type[q] == BD2) continue;
					coef.push_back(Tpl(p, q, 1.));
				}
				break;
			}
			d.setFromTriplets(coef.begin(), coef.end());
			mSol->a = mSol->a + d;
		}
#endif
		void makeDirchlet_v() {}

		void makeNeumann_p() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != BD1) continue;
				part->neumann[p] = 0.;
			}
		}

		void solvMat_p() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pres[p] = mSol->x[p];
				if (part->pres[p] < -1.e5) part->pres[p] = -1.e5;
				if (part->pres[p] > 1.e5) part->pres[p] = 1.e5;
			}
		}

		void solvMat_phi() {
			auto* const part = derived().part;
			mSol->ccBiCg_augment(part->type);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->phi[p] = mSol->x[p];
				if (part->phi[p] < -1.e5) part->phi[p] = -1.e5;
				if (part->phi[p] > 1.e5) part->phi[p] = 1.e5;
			}
		}

		void solvMat_v() {
			mSol->biCg_v();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				for (auto d = 0; d < D; d++) {
					part->vel2[p][d] = mSol->u[D*p + d];
				}
			}
		}

		void solvMat_t() {
			mSol->biCg();
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->temp[p] = mSol->x[p];
			}
		}

		R cfl() {
			R umax = 0.;
			const auto* const part = derived().part;
			for (unsigned p = 0; p < part->np; p++) {
				R tmp = part->vel1[p].norm();
				if (tmp > umax) umax = tmp;
			}
			para.umax = umax;
			return para.cfl * part->dp / umax;
		}

		void calBdNoSlip() {
			derived().part->bdNoSlip();
		}
		void bdSetZero() {
			derived().part->bdSetZero();
		}

		void calCell() {
			derived().part->updateCell();
		}

		void calPnd() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pnd[p] = part->cPnd(p);
				//part->pn[p] = part->cPn(p);
				//part->nbd[p] = part->cNbd(p);
			}
		}

		void calInvMat() {
			derived().part->updateInvMat();
		}

		void shift() {
			//shi.shiftPnd(part, para);
			//shi.shiftXu(derived().part, para);
			shi.shiftSpringIterate(derived().part, para);
			//shi.shiftNearest2d(derived().part, para);
			//shi.shiftPndLs(part, para);
			//shi.shiftLs(part, para);
			//shi.shiftCo(part, para);
			//shi.shiftPbf(part, para);
		}

		void makeFs() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) part->fs[p] = part->_isFs(p);
			//derived().part->updateTeam();
		}

		void pthOrderPresSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 1 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] == BD2 || part->isFs(p)) continue;
				part->pres[p] = part->func(part->pres, p);
			}
		}
		void pthOrderVelSpatialFilter() {
			static auto counter = 0;
			if (counter++ % 30 != 0) return;
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel2[p] = part->func(part->vel1, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID || part->isFs(p)) continue;
				part->vel1[p] = part->vel2[p];
			}
		}
		void surfCol() {
			auto* const part = derived().part;
			std::vector<vec> cor(part->np, vec(0.));
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID || !part->isFs(p)) continue;
				const auto c = part->cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < part->cell->linkList[key].size(); m++) {
						const auto q = part->cell->linkList[key][m];
						if (!part->isFs(q)) continue;
						if (q == p) continue;
						const auto dr = part->pos[q] - part->pos[p];
						const auto dv = part->vel1[q] - part->vel1[p];
						const auto dr1 = dr.mag();
						if (dr1 > part->r0) continue;
						const auto pro = dv*dr / dr1;
						if (abs(pro) > 0.8*para.umax) {
							cor[p] += 0.5* abs(pro)* dv.norm();
						}
					}
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				part->vel1[p] += cor[p];
				part->vel2[p] = part->vel1[p];
			}
		}
		void collision() {
			auto* const part = derived().part;
			std::vector<vec> cor(part->np, vec(0.));
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				const auto c = part->cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
						const auto q = part->cell->linkList[key][m];
						if (p == q) continue;
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.mag();
						if (dr1 < 0.5*part->dp) {
							const auto dv = part->vel1[q] - part->vel1[p];
							const auto tmp = dr*dv;
							if (tmp < 0.) {
								const auto coef = (0.5 / (dr1*dr1)*tmp*(1. + 0.2));
								cor[p] += coef* dr;
							}
						}
					}
				}
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel1[p] += cor[p];
			}
		}

		void damping() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID || !part->isFs(p)) continue;
				const auto n = part->vel1[p].norm();
				const auto gd = part->grad(part->vel1, p)* n;
				part->vel2[p] -= 0.01*part->dp* gd.mag()* n;
				part->vel1[p] = part->vel2[p];
			}
		}

		void calForVis() {
			auto* const part = derived().part;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(part->np); p++) {
				part->vort[p] = part->rot(part->vel1, p);
			}
		}

		void check() const {
			const auto* const part = derived().part;
			R velMax = std::numeric_limits<R>::min();
			R phiMax = std::numeric_limits<R>::min();
			R divMax = std::numeric_limits<R>::min();
			R divSum = R(0);
			unsigned idv = 0, idp = 0, idd = 0;
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) continue;
				const R vel = part->vel2[p].norm();
				const R phi = part->phi[p];
				const R div = part->div(part->vel2, p);
				if (vel > velMax) {
					velMax = vel;
					idv = p;
				}
				if (abs(phi) > abs(phiMax)) {
					phiMax = phi;
					idp = p;
				}
				if (abs(div) > abs(divMax)) {
					divMax = div;
					idd = p;
				}
				divSum += abs(div);
			}
			std::cout << " max vel: " << velMax << " --- id: " << idv << std::endl;
			std::cout << " max phi: " << phiMax << " --- id: " << idp << std::endl;
			std::cout << " max Div: " << divMax << " --- id: " << idd << std::endl;
			std::cout << " avg Div: " << divSum/part->np << std::endl;
		}

		void insertRand() {
			auto* const part = derived().part;
			R coef = 0.1;
			std::default_random_engine gen;
			std::normal_distribution<R> dis(0., 0.5);
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) continue;
				vec dr = coef* part->dp* vec(dis(gen), 0., dis(gen));
				part->pos[p] += dr;
			}
		}

	protected:
		int timeStep;

	};

}