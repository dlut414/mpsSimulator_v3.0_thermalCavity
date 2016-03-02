/*
*/
#pragma once
#include "Simulator.h"
#include "Particle_x_cst.h"

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class FractionalStep_DD : public Simulator <R, D, FractionalStep_DD<R, D, P>> {
		typedef mMath::Polynomial_A<R, D, P> PN;
		typedef mMath::Derivative_A<R, D, P> DR;
		typedef Eigen::Matrix<R, PN::value, 1> VecP;
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1>	Vec;
		typedef Eigen::Triplet<R>		Tpl;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
	public:
		FractionalStep_DD() {}
		~FractionalStep_DD() {}

		void step() {
			calInvMat();

			visTerm_i_q1();
			presTerm_i_q1();
			temperatureTerm_i_q1();

			syncPos();
			updateVelocity_q1();
			updatePressure_q1();
			updatePosition_s1();

			calInvMat();
			calForVis();
			check();

			calCell();
			calInvMat();
			shift();

			sync();
		}

		void visTerm_i_q1() {
			makeLhs_v_q1();
			makeRhs_v_q1();
			solvMat_v();
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID) part->vel2[p] += part->vel1[p];
			}
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

		void updateVelocity_q1() {
			const auto coefL = para.dt;
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
				if (part->type[p] == FLUID) part->pos[p] += para.dt * part->vel2[p];
			}
		}

		void updatePressure_q1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == FLUID || part->type[p] == BD1) part->pres[p] += part->phi[p];
			}
		}

		void init_() {
			part = new Particle_x_cst<R, D, P>();
			part->clean();
			*part << "Geo.in";
			part->init(para.k, para.beta);
			part->buildCell();
			part->makeBdc();
			part->b2b();
			part->b2norm();
			part->b2neumann();
			part->b2dirichlet();
			//part->updateTeam();
			part->init_x();
			sen = new Sensor<R, D, Particle_x_cst<R, D, P>>(part);
		}

	public:
		Particle_x_cst<R, D, P>* part;
		Sensor<R, D, Particle_x_cst<R, D, P>>* sen;

	private:
		__forceinline void makeLhs_v_q1() {
			coef.clear();
			for (unsigned p = 0; p < part->np; p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						coef.push_back(Tpl(D*p + d, D*p + d, 1.));
					}
					continue;
				}
				R pp = 0.;
				const R coefL = R(0.5)* para.dt* para.Pr;
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
						const auto pq = coefL* lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						for (auto d = 0; d < D; d++) {
							coef.push_back(Tpl(D*p + d, D*q + d, pq));
						}
					}
				}
				pp += R(1.);
				for (auto d = 0; d < D; d++) {
					coef.push_back(Tpl(D*p + d, D*p + d, pp));
				}
			}
			mSol->au.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeRhs_v_q1() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD1 || part->type[p] == BD2) {
					for (auto d = 0; d < D; d++) {
						mSol->rhs[D*p + d] = R(0.);
					}
					continue;
				}
				Vec g = Vec::Zero(); g[1] = para.Ra* part->temp[p];
				const Vec rhs = para.dt*( -part->grad(part->pres,p) + para.Pr* part->lap(part->vel1, p) + g );
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
				R pqsum = R(0.);
				R pp = R(0.);
				//std::vector<unsigned> used;
				MatPP* mm;
				if (IS(part->bdc[p], P_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
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

		__forceinline void makeRhs_p_q1() {
			const R coefL = R(1.) / para.dt;
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
					const auto cst = part->p_neumann.at(p)*part->w3(0.)* ((1. / part->varrho) * part->pn_lap_o) * (a);
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
				if (IS(part->bdc[p], T_DIRICHLET)) {
					coef.push_back(Tpl(p, p, 1.));
					continue;
				}
				R pqsum = 0.;
				R pp = 0.;
				MatPP* mm;
				if (IS(part->bdc[p], T_NEUMANN))	mm = &(part->invNeu.at(p));
				else								mm = &(part->invMat[p]);
				const R coefL = -R(0.5)* para.dt;
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
						const auto pq = coefL* lp.dot(a);
						pp -= pq;
						if (q == p) continue;
						coef.push_back(Tpl(p, q, pq));
						pqsum += abs(pq);
					}
				}
				pp += R(1.);
				coef.push_back(Tpl(p, p, pp));
				if (pqsum < para.eps) coef.push_back(Tpl(p, p, 1.));
			}
			mSol->a.setFromTriplets(coef.begin(), coef.end());
		}

		__forceinline void makeRhs_t_q1() {
			const R coefL = R(0.5)* para.dt;
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(part->np); p++) {
				if (part->type[p] == BD2) {
					mSol->b[p] = 0.;
					continue;
				}
				if (IS(part->bdc[p], T_DIRICHLET)) {
					mSol->b[p] = part->t_dirichlet.at(p);
					continue;
				}
				mSol->b[p] = part->temp[p] + coefL*part->lap(part->temp, p);
				if (IS(part->bdc[p], T_NEUMANN)) {
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