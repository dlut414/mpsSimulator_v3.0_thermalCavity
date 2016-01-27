/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_x_cst.h"

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class FractionalStep_x_cst : public Simulator <R,D,FractionalStep_x_cst<R,D,P>> {
		typedef mMath::Polynomial_A<R,D,P> PN;
		typedef mMath::Derivative_A<R,D,P> DR;
		typedef Eigen::Matrix<R,PN::value,1> VecP;
		typedef Eigen::Matrix<int,D,1>	iVec;
		typedef Eigen::Matrix<R,D,1>	Vec;
		typedef Eigen::Triplet<R>		Tpl;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
	public:
		FractionalStep_x_cst() {}
		~FractionalStep_x_cst() {}

		void step() {
			calInvMat();

			visTerm_i_q2r0();
			presTerm_i_q2();
			temperatureTerm_i_q1();

			syncPos();
			updateVelocity_q2r0();
			updatePosition_s2();

			calInvMat();
			calForVis();
			check();

			calCell();
			calInvMat();
			shift();

			sync();

			//if (timeStep % 100 == 0) profileOut_avgVel2();
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
			makeLhs_p();
			makeRhs_p_q2();
			solvMat_phi();
		}

		void presTerm_i_q1() {
			makeLhs_p();
			makeRhs_p_q1();
			solvMat_phi();
		}

		void temperatureTerm_i_q1() {
			makeLhs_t();
			makeRhs_t_q1();
			solvMat_t();
		}

		void updateVelocity_q2r1() {
			const auto coefL = (2.* para.dt) / (3.* para.rho);
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					part->pres[p] += part->phi[p] - miu* part->div(part->vel2, p);
				}
				else if (part->type[p] == BD2) part->pres[p] += part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->vel2[p] += -coefL* part->grad(part->phi, p);
			}
		}

		void updateVelocity_q1r0() {
			const auto coefL = para.dt / para.rho;
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID || part->type == BD1) {
					part->pres[p] = part->phi[p] - miu* part->div(part->vel2, p);
				}
				else if (part->type[p] == BD2) part->pres[p] = part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->vel2[p] += -coefL * part->grad(part->phi, p);
			}
		}

		void updateVelocity_q2r0() {
			const auto coefL = (2.* para.dt) / (3.* para.rho);
			const auto miu = para.rho * para.niu;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) {
					part->pres[p] = part->phi[p] - miu* part->div(part->vel2, p);
				}
				else if (part->type[p] == BD2) part->pres[p] = part->phi[p];
			}
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->vel2[p] += -coefL* part->grad(part->phi, p);
			}
		}

		void updatePosition_s1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->pos[p] += 0.5* para.dt * (part->vel1[p] + part->vel2[p]);
			}
		}

		void updatePosition_s2() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->pos[p] += 0.5* para.dt * (3.*part->vel1[p] - 1.* part->vel_m1[p]);
			}
		}

		void init_() {
			part = new Particle_x_cst<R,D,P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->b2b();
			part->b2norm();
			part->b2neumann();
			//part->updateTeam();
			part->init_x();
		}

	public:
		Particle_x_cst<R,D,P>* part;

	private:
		__forceinline void makeLhs_v_q2() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (int d = 0; d < D; d++) {
						coef.push_back(Tpl(D*p + d, D*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = mm * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = -para.niu* lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						for (auto d = 0; d < D; d++) {
							coef.push_back(Tpl(D*p + d, D*q + d, pq));
						}
					}
				}
				pp += 3. / (2. * para.dt);
				for (auto d = 0; d < D; d++) {
					coef.push_back(Tpl(D*p + d, D*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeLhs_v_q1() {
			coef.clear();
			for (auto p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						coef.push_back(Tpl(D*p + d, D*p + d, 1.));
					}
					continue;
				}
				auto pp = 0.;
				const auto& mm = part->invMat[p];
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.mag();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = mm * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = -para.niu* lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						for (auto d = 0; d < D; d++) {
							coef.push_back(Tpl(D*p + d, D*q + d, pq));
						}
					}
				}
				pp += 1. / para.dt;
				for (auto d = 0; d < D; d++) {
					coef.push_back(Tpl(D*p + d, D*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeRhs_v_q2r1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[D*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto gd = part->grad(part->pres, p);
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p])
					- (1. / para.rho)* gd
					+ para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[D*p + d] = rhs[d];
				}
			}
		}

		__forceinline void makeRhs_v_q1r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[D*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto rhs = (1. / para.dt)* part->vel1[p] + para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[D*p + d] = rhs[d];
				}
			}
		}

		__forceinline void makeRhs_v_q2r0() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[D*p + d] = part->vel1[p][d];
					}
					continue;
				}
				const auto rhs = 1. / (2.* para.dt)* (4.* part->vel1[p] - part->vel_m1[p]) + para.g;
				for (auto d = 0; d < D; d++) {
					mSol->rhs[D*p + d] = rhs[d];
				}
			}
		}


		__forceinline void makeLhs_p() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				if (part->isFs(p)) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				auto pqsum = 0.;
				auto pp = 0.;
				//std::vector<unsigned> used;
				MatPP* mm;
				if (part->type[p] != BD1)	mm = &(part->invMat[p]);
				else						mm = &(part->invNeu[p]);
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				//if (abs(mm.determinant()) > part->eps_mat) {}
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					//for (unsigned us = 0; us < used.size(); us++) { if (key == used[us]) std::cout << "used!!!!!!!!!!" << std::endl; }
					//used.push_back(key);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
#if NOBD2
						if (part->type[q] == BD2) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = (*mm) * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
						pqsum += abs(pq);
					}
				}
				coef.push_back(Tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, 1.));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeRhs_p_q2() {
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
				if (part->type[p] == BD1) {
					VecP inner = VecP::Zero();
					inner.block<D, 1>(0, 0) = part->bdnorm.at(p);
					const auto a = part->invMat[p]* inner;
					const auto cst = part->neumann.at(p)*part->w3(0.)* ((1./part->varrho) * part->pn_lap_o) * (a);
					mSol->b[p] -= cst;
				}
			}
		}

		__forceinline void makeRhs_p_q1() {
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
				if (part->type[p] == BD1) {
					VecP inner = VecP::Zero();
					inner.block<D, 1>(0, 0) = part->bdnorm.at(p);
					const auto a = part->invMat[p] * inner;
					const auto cst = part->neumann.at(p)*part->w3(0.)* ((1. / part->varrho) * part->pn_lap_o) * (a);
					mSol->b[p] -= cst;
				}
			}
		}

		__forceinline void makeLhs_t() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD2) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				if (part->bdc & T_DIRICHLET) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				R pqsum = 0.;
				R pp = 0.;
				MatPP* mm;
				if (part->bdc & T_NEUMANN)	mm = &(part->invMat[p]);
				else						mm = &(part->invNeu[p]);
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (part->bdOpt(p, q)) continue;
#endif
#if NOBD2
						if (part->type[q] == BD2) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						const auto w = part->w3(dr1);
						VecP npq;
						part->poly(dr, npq);
						const auto a = (*mm) * (w* npq);
						const auto& lp = part->pn_lap_o;
						const auto pq = lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
						pqsum += abs(pq);
					}
				}
				pp += -(1. / para.dt);
				coef.push_back(Tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, 1.));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeRhs_t_q1() {
			const auto coefL = -(1. / para.dt);
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (part->bdc & T_DIRICHLET) {
					mSol->b[p] = part->t_dirichlet.at(p);
					continue;
				}
				mSol->b[p] = coefL * part->temp[p];
				if (part->bdc & T_NEUMANN) {
					VecP inner = VecP::Zero();
					inner.block<D, 1>(0, 0) = part->bdnorm.at(p);
					const auto a = part->invMat[p] * inner;
					const auto cst = part->t_neumann.at(p)*part->w3(0.)* ((1. / part->varrho) * part->pn_lap_o) * (a);
					mSol->b[p] -= cst;
				}
			}
		}


		__forceinline void sync() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->vel_m1[p] = part->vel1[p];
				part->vel1[p] = part->vel2[p];
			}
		}
		__forceinline void syncPos() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				part->pos_m1[p] = part->pos[p];
			}
		}

	private:
		std::vector<Tpl> coef;
	};

}