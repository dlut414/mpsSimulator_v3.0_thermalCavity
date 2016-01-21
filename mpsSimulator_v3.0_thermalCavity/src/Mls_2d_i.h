/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_2d_i.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class Mls_2d_i : public Simulator < real, dim, Mls_2d_i<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Matrix<real, 5, 1> vecp;
		typedef Eigen::Triplet<real> tpl;
	public:
		Mls_2d_i() {}
		~Mls_2d_i() {}

		void step() {
			calInvMat();

			visTerm_i_q2r0();
			calPnd();
			presTerm_i_q2();

			convect_q2r0s2();

			calInvMat();
			calForVis();
			check();

			calCell(); //update cells
			calInvMat();
			shift();

			sync();

			if (timeStep % 100 == 0) profileOut_avgVel2();
			//calInvMat(); //for sensor
			//profileOut();
			//sensorOut();
			//pthOrderVelSpatialFilter();
		}

		void visTerm_e() {
			/*Euler*/
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] = part->vel1[p] + para.dt * (para.g + para.niu * part->lap(part->vel1, p));
			}
		}

		void visTerm_i_q2r1() {
			makeLhs_v_q2();
			makeRhs_v_q2r1();
			solvMat_v();
		}

		void visTerm_i_q1r0() {
			makeLhs_v_q1();
			makeRhs_v_q1r0();
			solvMat_v();
		}

		void visTerm_i_q2r0() {
			makeLhs_v_q2();
			makeRhs_v_q2r0();
			solvMat_v();
		}

		void presTerm_e() {}

		void presTerm_i_q2() {
			//makeDirichlet_p_op();
			makeLhs_p();
			makeRhs_p_q2();
			//makeDirichlet_p_avg();
			solvMat_phi();
		}

		void presTerm_i_q1() {
			//makeDirichlet_p_op();
			makeLhs_p();
			makeRhs_p_q1();
			//makeDirichlet_p_avg();
			solvMat_phi();
		}

		void convect_q2r1s2() {
			const auto coefL = (2.* para.dt) / (3.* para.rho);
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != BD2) {
					part->pres[p] += part->phi[p] - miu* part->div(part->vel2, p);
				}
				else part->pres[p] += part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -coefL* part->grad(part->phi, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (3.*part->vel1[p] - 1.* part->vel_m1[p]);
			}
		}

		void convect_q1r0s1() {
			const auto coefL = para.dt / para.rho;
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != BD2) {
					part->pres[p] = part->phi[p] - miu* part->div(part->vel2, p);
				}
				else part->pres[p] = part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -coefL * part->grad(part->phi, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
				part->vel1[p] = part->vel2[p];
			}
		}

		void convect_q2r0s1() {
			const auto coefL = (2.* para.dt) / (3.* para.rho);
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != BD2) {
					part->pres[p] = part->phi[p] - miu* part->div(part->vel2, p);
				}
				else part->pres[p] = part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -coefL* part->grad(part->phi, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
			}
		}

		void convect_q2r0s2() {
			const auto coefL = (2.* para.dt) / (3.* para.rho);
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != BD2) {
					part->pres[p] = part->phi[p] - miu* part->div(part->vel2, p);
				}
				else part->pres[p] = part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->vel2[p] += -coefL* part->grad(part->phi, p);
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) continue;
				part->pos[p] += 0.5* para.dt * (3.*part->vel1[p] - 1.* part->vel_m1[p]);
			}
		}

		void init_() {
			part = new Particle_2d_i<real, dim>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->b2b();
			part->b2norm();
			//part->updateTeam();
			part->init2d_x();
		}

	public:
		Particle_2d_i<real, dim>*  part;

	private:
		void makeLhs_v_q2() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) {
					for (int d = 0; d < dim; d++) {
						coef.push_back(tpl(dim*p + d, dim*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto mm = part->invMat[p];
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
#if BD_OPT
								if (part->bdOpt(p, q)) continue;
#endif
								const auto	dr = part->pos[q] - part->pos[p];
								const auto	dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								const auto w = part->w3(dr1);
								const auto npq = part->poly(dr);
								const auto a = mm * (w* npq);
								const auto lp = part->poly_lap_0;
								const auto pq = -para.niu* lp.dot(a);
								pp -= pq;
								if (q == p) continue;
								for (int d = 0; d < dim; d++) {
									coef.push_back(tpl(dim*p + d, dim*q + d, pq));
								}
							}
						}
					}
				}
				pp += 3. / (2. * para.dt);
				for (int d = 0; d < dim; d++) {
					coef.push_back(tpl(dim*p + d, dim*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		void makeLhs_v_q1() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] != FLUID) {
					for (int d = 0; d < dim; d++) {
						coef.push_back(tpl(dim*p + d, dim*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto mm = part->invMat[p];
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
#if BD_OPT
								if (part->bdOpt(p, q)) continue;
#endif
								const auto	dr = part->pos[q] - part->pos[p];
								const auto	dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								if (q == p) continue;
								const auto w = part->w3(dr1);
								const auto npq = part->poly(dr);
								const auto a = mm * (w* npq);
								const auto lp = part->poly_lap_0;
								const auto pq = -para.niu* lp.dot(a);
								pp -= pq;
								for (int d = 0; d < dim; d++) {
									coef.push_back(tpl(dim*p + d, dim*q + d, pq));
								}
							}
						}
					}
				}
				pp += 1. / para.dt;
				for (int d = 0; d < dim; d++) {
					coef.push_back(tpl(dim*p + d, dim*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}


		void makeRhs_v_q2r1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) {
					switch (dim) {
					case TWOD:
						mSol->rhs[2 * p] = part->vel1[p].x;
						mSol->rhs[2 * p + 1] = part->vel1[p].z;
						break;
					case THREED:
						mSol->rhs[3 * p] = part->vel1[p].x;
						mSol->rhs[3 * p + 1] = part->vel1[p].y;
						mSol->rhs[3 * p + 2] = part->vel1[p].z;
						break;
					default:
						break;
					}
					continue;
				}
				const auto gd = part->grad(part->pres, p);
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p])
					- (1. / para.rho)* gd
					+ para.g;
				switch (dim) {
				case TWOD:
					mSol->rhs[2 * p] = rhs.x;
					mSol->rhs[2 * p + 1] = rhs.z;
					break;
				case THREED:
					mSol->rhs[3 * p] = rhs.x;
					mSol->rhs[3 * p + 1] = rhs.y;
					mSol->rhs[3 * p + 2] = rhs.z;
					break;
				default:
					break;
				}
			}
		}

		void makeRhs_v_q1r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) {
					switch (dim) {
					case TWOD:
						mSol->rhs[2 * p] = part->vel1[p].x;
						mSol->rhs[2 * p + 1] = part->vel1[p].z;
						break;
					case THREED:
						mSol->rhs[3 * p] = part->vel1[p].x;
						mSol->rhs[3 * p + 1] = part->vel1[p].y;
						mSol->rhs[3 * p + 2] = part->vel1[p].z;
						break;
					default:
						break;
					}
					continue;
				}
				const auto rhs = (1. / para.dt)* part->vel1[p] + para.g;
				switch (dim) {
				case TWOD:
					mSol->rhs[2 * p] = rhs.x;
					mSol->rhs[2 * p + 1] = rhs.z;
					break;
				case THREED:
					mSol->rhs[3 * p] = rhs.x;
					mSol->rhs[3 * p + 1] = rhs.y;
					mSol->rhs[3 * p + 2] = rhs.z;
					break;
				default:
					break;
				}
			}
		}

		void makeRhs_v_q2r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] != FLUID) {
					switch (dim) {
					case TWOD:
						mSol->rhs[2 * p] = part->vel1[p].x;
						mSol->rhs[2 * p + 1] = part->vel1[p].z;
						break;
					case THREED:
						mSol->rhs[3 * p] = part->vel1[p].x;
						mSol->rhs[3 * p + 1] = part->vel1[p].y;
						mSol->rhs[3 * p + 2] = part->vel1[p].z;
						break;
					default:
						break;
					}
					continue;
				}
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p]) + para.g;
				switch (dim) {
				case TWOD:
					mSol->rhs[2 * p] = rhs.x;
					mSol->rhs[2 * p + 1] = rhs.z;
					break;
				case THREED:
					mSol->rhs[3 * p] = rhs.x;
					mSol->rhs[3 * p + 1] = rhs.y;
					mSol->rhs[3 * p + 2] = rhs.z;
					break;
				default:
					break;
				}
			}
		}


		void makeLhs_p() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(tpl(p, p, 1.));
					coef.push_back(tpl(p, part->bbMap.at(p), -1.));
					continue;
				}
				if (part->isFs(p)) {
					coef.push_back(tpl(p, p, 1.));
					continue;
				}
				auto pqsum = 0.;
				auto pp = 0.;
				//std::vector<unsigned> used;
				const auto mm = part->invMat[p];
				const auto c = part->cell->iCoord(part->pos[p]);
				//if (abs(mm.determinant()) > part->eps_mat) {
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							//for (unsigned us = 0; us < used.size(); us++) { if (key == used[us]) std::cout << "used!!!!!!!!!!" << std::endl; }
							//used.push_back(key);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const auto q = part->cell->linkList[key][m];
#if BD_OPT
								if (part->bdOpt(p, q)) continue;
#endif
								const auto	dr = part->pos[q] - part->pos[p];
								const auto	dr1 = dr.mag();
								if (dr1 > part->r0) continue;
								const auto w = part->w3(dr1);
								const auto npq = part->poly(dr);
								const auto a = mm * (w* npq);
								const auto lp = part->poly_lap_0;
								const auto pq = lp.dot(a);
								pp -= pq;
								if (q == p) continue;
								coef.push_back(tpl(p, q, pq));
								pqsum += abs(pq);
							}
						}
					}
				}
				coef.push_back(tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(tpl(p, p, 1.));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeRhs_p_q2() {
			const auto coefL = (3.* para.rho) / (2.* para.dt);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = 0.;
					continue;
				}
				mSol->b[p] = coefL * part->div(part->vel2, p);
				//mSol->b[p] = coefL * ( part->div(part->vel2, p) + (0.3/para.dt)*(part->n0-part->pnd[p])/part->n0 );
			}
		}

		void makeRhs_p_q1() {
			const auto coefL = para.rho / para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->isFs(p)) {
					mSol->b[p] = 0.;
					continue;
				}
				mSol->b[p] = coefL * part->div(part->vel2, p);
			}
		}

		
		void sync() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->vel_m1[p] = part->vel1[p];
				part->vel1[p] = part->vel2[p];
			}
		}

	private:
		std::vector<tpl> coef;
	};

}