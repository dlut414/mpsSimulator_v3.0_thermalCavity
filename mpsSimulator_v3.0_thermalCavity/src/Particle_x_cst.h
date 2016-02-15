/*
*/
#pragma once
#include "Header.h"
#include "Particle.h"
#include "Polynomial.h"
#include "Derivative.h"

#define UPWIND_VEL	0
#define NOBD2		1

namespace SIM {

	template <typename R, unsigned D, unsigned P>
	class Particle_x_cst : public Particle<R,D,Particle_x_cst<R,D,P>> {
		typedef mMath::Polynomial_A<R, D, P> PN;
		typedef mMath::Derivative_A<R, D, P> DR;
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1>		Vec;
		typedef Eigen::Matrix<R, D, D>		Mat;
		typedef Eigen::Matrix<R, PN::value, 1>	VecP;
		typedef Eigen::Matrix<R, PN::value, D>	MatPD;
		typedef Eigen::Matrix<R, PN::value, PN::value> MatPP;
	public:
		Particle_x_cst() : Particle() {}
		~Particle_x_cst() {}
		
		__forceinline void poly(const Vec& in, VecP& out) const { PN::Gen(varrho, in.data(), out.data()); }

		const Vec grad(const std::vector<R>& phi, const unsigned& p) const {
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += w * (phi[q] - phi[p]) * npq;
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		const Mat grad(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w* npq)* (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_p_o*a);
		}

		template <typename T, typename U, typename V>
		const T grad(const U& phi, const V& p) const {
		}

		const R div(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w* npq)* (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			auto ret = static_cast<R>(0);
			for (auto d = 0; d < D; d++) {
				ret += pn_p_o.block<1, PN::value>(d, 0) * a.block<PN::value, 1>(0, d);
			}
			return ret;
		}

		const R lap(const std::vector<R>& phi, const unsigned& p) const {
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP	npq;
					poly(dr, npq);
					vv += w * (phi[q] - phi[p])* npq;
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_lap_o*a);
		}

		const Vec lap(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w * npq) * (u[q] - [p]);
				}
			}
			const auto a = invMat[p] * vv;
			return (pn_lap_o*a).transpose();
		}

		const R rot(const std::vector<Vec>& u, const unsigned& p) const {
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					vv += (w * npq) * (u[q] - u[p]).transpose();
				}
			}
			const auto a = invMat[p] * vv;
			const auto der = pn_p_o*a;
			switch (D) {
			case 1:
				return R(0);
				break;
			case 2:
				return R(der(0, 1) - der(1, 0));
				break;
			case 3:
				return R(0);
				break;
			default:
				return R(0);
			}
		}

		const R func(const std::vector<R>& phi, const unsigned& p) const {
			return phi[p];
		}

		const Vec func(const std::vector<Vec>& u, const unsigned& p) const {
			return u[p];
		}

		const R func(const std::vector<R>& phi, const Vec& p) const {
			auto rid = 0;
			auto isNear = 0;
			auto rr = std::numeric_limits<R>::max();
			auto c = cell->iCoord(p);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
					const auto dr = pos[q] - p;
					const auto dr1 = dr.norm();
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					if (dr1 > r0) continue;
					else {
						isNear = 1;
						if (dr1 < rr) {
							rr = dr1;
							rid = q;
						}
					}
				}
			}
			if (!isNear) return R(0);
			else return phi[rid] + (p - pos[rid]).transpose()*grad(phi, rid);
		}

		const Vec func(const std::vector<Vec>& phi, const Vec& p) const {
			auto rid = 0;
			auto isNear = 0;
			auto rr = std::numeric_limits<R>::max();
			auto c = cell->iCoord(p);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
					const auto dr = pos[q] - p;
					const auto dr1 = dr.norm();
#if NOBD2
					if (type[q] == BD2) continue;
#endif
					if (dr1 > r0) continue;
					else {
						isNear = 1;
						if (dr1 < rr) {
							rr = dr1;
							rid = q;
						}
					}
				}
			}
			if (!isNear) return Vec::Zero();
			else {
				const auto dpt = (p - pos[rid]).transpose();
				return phi[rid] + (dpt*grad(phi, rid)).transpose();
			}
		}

		const Vec func_mafl(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (auto i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = phi[p];
			for (auto fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec::Zero();
				auto ww = static_cast<R>(0);
				auto c = cell->iCoord(pLocal[fp]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
						if (type[q] == BD2) continue;
#if NOBD2
						if (type[q] == BD2) continue;
#endif
#if BD_OPT
						if (bdOpt(q)) continue;
#endif
						const auto dr1 = (pos[q] - pLocal[fp]).norm();
						const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).norm();
						const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).norm();
						if (dr1 > re) continue;
						if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
						const auto w = w1(dr1);
						ww += w;
						ret[fp] += w * phi[q];
					}
				}
				if (abs(ww) < eps) ww = 1.;
				ret[fp] = ret[fp] / ww;
			}
			return ret[3] - (dmove.norm() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
		}

		const Vec func_mafl_mmt(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			Vec pLocal[6];

			for (auto i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			Vec ret[6];
			ret[3] = phi[p];
			for (auto fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = Vec::Zero();
				auto ww = static_cast<R>(0);
				auto c = cell->iCoord(pLocal[fp]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
						if (type[q] == BD2) continue;
#if NOBD2
						if (type[q] == BD2) continue;
#endif
#if BD_OPT
						if (bdOpt(q)) continue;
#endif
						const auto dr1 = (pos[q] - pLocal[fp]).norm();
						const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).norm();
						const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).norm();
						if (dr1 > re) continue;
						if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
						const auto w = w1(dr1);
						ww += w;
						ret[fp] += w * phi[q];
					}
				}
				if (abs(ww) < eps) continue;
				ret[fp] = ret[fp] / ww;
			}
			auto ret_mmt = ret[3] - (dmove.norm() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
			auto ret_min = ret[1];
			auto ret_max = ret[1];
			for (int fp = 1; fp <= 4; fp++) {
				if (ret[fp].squaredNorm() < ret_min.squaredNorm()) ret_min = ret[fp];
				if (ret[fp].squaredNorm() > ret_max.squaredNorm()) ret_max = ret[fp];
			}
			if (ret_mmt.squaredNorm() < ret_min.squaredNorm()) ret_mmt = ret_min* ret_mmt.norm();
			if (ret_mmt.squaredNorm() > ret_max.squaredNorm()) ret_mmt = ret_max* ret_mmt.norm();
			return ret_mmt;
		}

		const Vec func_lsA(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			auto mm = MatPP::Zero();
			auto vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if NOBD2
					if (type[q] == BD2) continue;
#endif
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]).transpose();
				}
			}
			const auto inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					inv = MatPP::Zero();
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
			const Mat hes[D];
			for (auto d = 0; d < D; d++) {
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes[d](i, j) = mgd.block<mMath::H<D,2>,1>(0, d);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < D; j++) {
						hes[d](i, j) = hes[d](j, i);
					}
				}
			}
			const auto dp = p_new - pos[p];
			const auto dpt = dp.transpose();
			auto ret = phi[p];
			for (auto d = 0; d < D; d++) ret[d] += (dpt*gd).transpose() + 0.5*dpt*hes[d]*dp;
			return ret;
		}

		const Vec func_lsA_upwind(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(phi[p].norm());
#else
			const auto up = dp.normalized();
#endif
			MatPP mm = MatPP::Zero();
			MatPD vv = MatPD::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if NOBD2
					if (type[q] == BD2) continue;
#endif
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					if (dr.dot(up) < 0) continue;
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* npq)* (phi[q] - phi[p]).transpose();
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					//inv = MatPP::Zero();
					inv = invMat[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
			Mat hes[D];
			int counter = 0;
			for (auto d = 0; d < D; d++) {
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes[d](i, j) = mgd(counter++, d);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < i; j++) {
						hes[d](i, j) = hes[d](j, i);
					}
				}
			}
			const auto dpt = dp.transpose();
			auto ret = phi[p];
			for (auto d = 0; d < D; d++) ret[d] += (dpt*gd)[d] + 0.5*dpt * hes[d] * dp;
			return ret;
		}

		const R func_lsA_upwind(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(phi[p].norm());
#else
			const auto up = dp.normalized();
#endif
			MatPP mm = MatPP::Zero();
			VecP vv = VecP::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (auto m = 0; m < cell->linkList[key].size(); m++) {
					const auto q = cell->linkList[key][m];
#if NOBD2
					if (type[q] == BD2) continue;
#endif
#if BD_OPT
					if (bdOpt(p, q)) continue;
#endif
					const auto dr = pos[q] - pos[p];
					if (dr.dot(up) < 0) continue;
					const auto dr1 = dr.norm();
					if (dr1 > r0) continue;
					const auto w = w3(dr1);
					VecP npq;
					poly(dr, npq);
					mm += (w* npq)* npq.transpose();
					vv += (w* (phi[q] - phi[p]))* npq;
				}
			}
			MatPP inv = MatPP::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<2, 2>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					//inv = MatPP::Zero();
					inv = invMat[p];
				}
				else inv.block<2, 2>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto gd = pn_p_o * a;
			const auto mgd = pn_pp_o * a;
			Mat hes;
			int counter = 0;
			for (auto i = 0; i < D; i++) {
				for (auto j = i; j < D; j++) {
					hes(i, j) = mgd(counter++);
				}
			}
			for (auto i = 0; i < D; i++) {
				for (auto j = 0; j < i; j++) {
					hes(i, j) = hes(j, i);
				}
			}
			const auto dpt = dp.transpose();
			R ret = phi[p];
			ret = ret + dpt* gd + 0.5* dpt * hes* dp;
			return ret;
		}

		//const Vec func_lsB(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}

		//const Vec func_lsB_upwind(const std::vector<Vec>& u, const unsigned& p, const Vec& p_new) const {}

		template <int StencilsX = 1, int StencilsY = 3, int Stencils = StencilsX*StencilsY, int Dimension = D>	struct interpolateWENO_B_ {
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 1> {
			template <typename R> static const R Gen(const std::vector<R>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {}
			template <typename Vec> static const Vec Gen(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {}
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 2> {
			template <typename R> static const R Gen(const std::vector<R>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {
				const auto dp = p_new - part->pos[p];
				if (dp.norm() < part->eps) return phi[p];
				const auto up = dp.normalized();
				const auto alpha = 2.* M_PI / StencilsX;
				Vec dir[StencilsX];
				Vec ctr[Stencils];
				for (auto i = 0; i < StencilsX; i++) {
					const auto theta = i* alpha;
					const auto ct = cos(theta);
					const auto st = sin(theta);
					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
				}
				for (auto j = 0; j < StencilsY; j++) {
					//const auto dis = part->r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
					const auto dis = part->r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
					for (auto i = 0; i < StencilsX; i++) {
						const auto stcId = i* StencilsY + j;
						ctr[stcId] = part->pos[p] + dis*dir[i];
					}
				}
				MatPP mm[Stencils];
				VecP vv[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					mm[i] = MatPP::Zero();
					vv[i] = VecP::Zero();
				}
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						for (auto stcX = 0; stcX < StencilsX; stcX++) {
							for (auto stcY = 0; stcY < StencilsY; stcY++) {
								const auto stcId = stcX* StencilsY + stcY;
								const auto dis = (part->pos[q] - ctr[stcId]).norm();
								const auto w = part->w3(dis);
								VecP npq;
								part->poly(dr, npq);
								mm[stcId] += (w* npq)* npq.transpose();
								vv[stcId] += (w* npq)* (phi[q] - phi[p]);
							}
						}
					}
				}
				R oscillationIndicator[Stencils];
				R stencilWeight[Stencils];
				R stencilWeightNorm[Stencils];
				VecP polyCoef[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					MatPP inv = MatPP::Zero();
					if (abs(mm[i].determinant()) < part->eps_mat) {
						return phi[p];
						auto mm_ = mm[i].block<2, 2>(0, 0);
						if (abs(mm_.determinant()) < part->eps_mat) {
							inv = MatPP::Zero();
						}
						else inv.block<2, 2>(0, 0) = mm_.inverse();
					}
					else inv = mm[i].inverse();
					polyCoef[i] = inv * vv[i];
					oscillationIndicator[i] = R(0.);
				}

				for (auto i = 0; i < Stencils; i++) {
					const auto dp = part->dp;
					const auto A = polyCoef[i][0] * polyCoef[i][0];
					const auto B = polyCoef[i][1] * polyCoef[i][1];
					const auto C = polyCoef[i][2] * polyCoef[i][2];
					const auto D = polyCoef[i][3] * polyCoef[i][3];
					const auto E = polyCoef[i][4] * polyCoef[i][4];
					const auto beta1 = (dp*dp*dp)*(A + B) + (dp*dp*dp*dp*dp / 6.)*(2.*C + D + 2.*E);
					const auto beta2 = (dp*dp*dp*dp*dp)*(4.*C + D + 4.*E);
					//oscillationIndicator[i] = 4.*beta1 - (1. / 3.)*beta2;
					oscillationIndicator[i] = beta1 + beta2;
				}
				const R epsilon = 1.e-6;
				const int magnifier = 5;
				for (auto i = 0; i < Stencils; i++) {
					stencilWeight[i] = 1. / pow(epsilon + oscillationIndicator[i], magnifier);
				}
				R stencilWeightSum = R(0.);
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightSum += stencilWeight[i];
				}
				for (auto i = 0; i < Stencils; i++) {
					stencilWeightNorm[i] = stencilWeight[i] / stencilWeightSum;
				}
				VecP combinedCoef = VecP::Zero();
				for (auto i = 0; i < Stencils; i++) {
					combinedCoef += stencilWeightNorm[i] * polyCoef[i];
				}
				const auto gd = part->pn_p_o * combinedCoef;
				const auto mgd = part->pn_pp_o * combinedCoef;
				Mat hes;
				int counter = 0;
				for (auto i = 0; i < D; i++) {
					for (auto j = i; j < D; j++) {
						hes(i, j) = mgd(counter++);
					}
				}
				for (auto i = 0; i < D; i++) {
					for (auto j = 0; j < i; j++) {
						hes(i, j) = hes(j, i);
					}
				}
				const auto dpt = dp.transpose();
				auto ret = phi[p];
				ret = ret + (dpt*gd) + 0.5*dpt * hes * dp;
				return ret;
			}
			template <typename Vec> static const Vec Gen(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {
				const auto dp = p_new - part->pos[p];
				if (dp.norm() < part->eps) return phi[p];
				const auto up = dp.normalized();
				const auto alpha = 2.* M_PI / StencilsX;
				Vec dir[StencilsX];
				Vec ctr[Stencils];
				for (auto i = 0; i < StencilsX; i++) {
					const auto theta = i* alpha;
					const auto ct = cos(theta);
					const auto st = sin(theta);
					dir[i] << ct*up[0] + st*up[1], ct*up[1] - st*up[0];
				}
				for (auto j = 0; j < StencilsY; j++) {
					//const auto dis = part->r0* ( R(1.) - R(2.)*(j + 1) / (1 + StencilsY) );
					const auto dis = part->r0* (R(1.) - R(1.)*(j + 1) / (StencilsY));
					for (auto i = 0; i < StencilsX; i++) {
						const auto stcId = i* StencilsY + j;
						ctr[stcId] = part->pos[p] + dis*dir[i];
					}
				}
				MatPP mm[Stencils];
				MatPD vv[Stencils];
				for (auto i = 0; i < Stencils; i++) {
					mm[i] = MatPP::Zero();
					vv[i] = MatPD::Zero();
				}
				const auto& cell = part->cell;
				const auto c = cell->iCoord(part->pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if BD_OPT
						if (bdOpt(p, q)) continue;
#endif
						const auto dr = part->pos[q] - part->pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > part->r0) continue;
						for (auto stcX = 0; stcX < StencilsX; stcX++) {
							for (auto stcY = 0; stcY < StencilsY; stcY++) {
								const auto stcId = stcX* StencilsY + stcY;
								const auto dis = (part->pos[q] - ctr[stcId]).norm();
								const auto w = part->w3(dis);
								VecP npq;
								part->poly(dr, npq);
								mm[stcId] += (w* npq)* npq.transpose();
								vv[stcId] += (w* npq)* (phi[q] - phi[p]).transpose();
							}
						}
					}
				}
				R oscillationIndicator[D][Stencils];
				R stencilWeight[D][Stencils];
				R stencilWeightNorm[D][Stencils];
				VecP polyCoef[D][Stencils];
				for (auto d = 0; d < D; d++) {
					for (auto i = 0; i < Stencils; i++) {
						MatPP inv = MatPP::Zero();
						if (abs(mm[i].determinant()) < part->eps_mat) {
							return phi[p];
							auto mm_ = mm[i].block<2, 2>(0, 0);
							if (abs(mm_.determinant()) < part->eps_mat) {
								inv = MatPP::Zero();
							}
							else inv.block<2, 2>(0, 0) = mm_.inverse();
						}
						else inv = mm[i].inverse();
						polyCoef[d][i] = (inv * vv[i]).block<PN::value, 1>(0, d);
						oscillationIndicator[d][i] = R(0.);
					}
				}

				for (auto d = 0; d < D; d++) {
					for (auto i = 0; i < Stencils; i++) {
						const auto dp = part->dp;
						const auto A = polyCoef[d][i][0] * polyCoef[d][i][0];
						const auto B = polyCoef[d][i][1] * polyCoef[d][i][1];
						const auto C = polyCoef[d][i][2] * polyCoef[d][i][2];
						const auto D = polyCoef[d][i][3] * polyCoef[d][i][3];
						const auto E = polyCoef[d][i][4] * polyCoef[d][i][4];
						const auto beta1 = (dp*dp*dp)*(A + B) + (dp*dp*dp*dp*dp / 6.)*(2.*C + D + 2.*E);
						const auto beta2 = (dp*dp*dp*dp*dp)*(4.*C + D + 4.*E);
						//oscillationIndicator[d][i] = 4.*beta1 - (1. / 3.)*beta2;
						oscillationIndicator[d][i] = beta1 + beta2;
					}
				}
				const R epsilon = 1.e-6;
				const int magnifier = 5;
				for (auto d = 0; d < D; d++) {
					for (auto i = 0; i < Stencils; i++) {
						stencilWeight[d][i] = 1. / pow(epsilon + oscillationIndicator[d][i], magnifier);
					}
					R stencilWeightSum = R(0.);
					for (auto i = 0; i < Stencils; i++) {
						stencilWeightSum += stencilWeight[d][i];
					}
					for (auto i = 0; i < Stencils; i++) {
						stencilWeightNorm[d][i] = stencilWeight[d][i] / stencilWeightSum;
					}
				}
				auto ret = phi[p];
				for (auto d = 0; d < D; d++) {
					VecP combinedCoef = VecP::Zero();
					for (auto i = 0; i < Stencils; i++) {
						combinedCoef += stencilWeightNorm[d][i] * polyCoef[d][i];
					}
					const auto gd = part->pn_p_o * combinedCoef;
					const auto mgd = part->pn_pp_o * combinedCoef;
					Mat hes;
					int counter = 0;
					for (auto i = 0; i < D; i++) {
						for (auto j = i; j < D; j++) {
							hes(i, j) = mgd(counter++);
						}
					}
					for (auto i = 0; i < D; i++) {
						for (auto j = 0; j < i; j++) {
							hes(i, j) = hes(j, i);
						}
					}
					const auto dpt = dp.transpose();
					ret[d] = ret[d] + (dpt*gd) + 0.5*dpt * hes * dp;
				}
				return ret;
			}
		};
		template <int StencilsX, int StencilsY, int Stencils>		struct interpolateWENO_B_<StencilsX, StencilsY, Stencils, 3> {
			template <typename R> static const R Gen(const std::vector<R>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {}
			template <typename Vec> static const Vec Gen(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new, Particle_x_cst<R, D, P>* const part) {}
		};

		__forceinline const R interpolateWENO(const std::vector<R>& phi, const unsigned& p, const Vec& p_new) {
			return interpolateWENO_B_<>::Gen(phi, p, p_new, this);
		}
		__forceinline const Vec interpolateWENO(const std::vector<Vec>& phi, const unsigned& p, const Vec& p_new) {
			return interpolateWENO_B_<>::Gen(phi, p, p_new, this);
		}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				MatPP mm = MatPP::Zero();
				const auto c = cell->iCoord(pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (auto m = 0; m < cell->linkList[key].size(); m++) {
						const auto q = cell->linkList[key][m];
#if NOBD2
						if (type[q] == BD2) continue;
#endif
#if BD_OPT
						if (bdOpt(p, q)) continue;
#endif
						const auto dr = pos[q] - pos[p];
						const auto dr1 = dr.norm();
						if (dr1 > r0) continue;
						const auto w = w3(dr1);
						VecP npq;
						poly(dr, npq);
						mm += (w* npq) * npq.transpose();
					}
				}

				if (type[p] == BD1) {
					const auto& inner = bdnorm.at(p);
					MatPP nn = MatPP::Zero();
					nn.block<D, D>(0, 0) = w3(0)* inner * inner.transpose();
					MatPP mpn = mm + nn;
					auto& inv = invNeu.at(p);
					inv = MatPP::Zero();
					if (abs(mpn.determinant()) < eps_mat) {
#if DEBUG
						std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mpn.determinant() << std::endl;
#endif
						auto mpn_ = mpn.block<2, 2>(0, 0);
						if (abs(mpn_.determinant()) < eps_mat) inv = MatPP::Zero();
						else inv.block<2, 2>(0, 0) = mpn_.inverse();
					}
					else inv = mpn.inverse();
				}

				auto& invRef = invMat.at(p);
				invRef = MatPP::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << " ID: " << p << " --- " << " Determinant defficiency: " << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<2, 2>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) invRef = MatPP::Zero();
					else invRef.block<2, 2>(0, 0) = mm_.inverse();
				}
				else invRef = mm.inverse();
			}
		}

		template <unsigned D_ = D>
		void init_x() {}

		template <>
		void init_x<1>() {
			varrho = 1./(1.*dp);
			Vec zero = Vec::Zero();
			DR::Gen<1>(varrho, zero.data(), pn_p_o.data());
			DR::Gen<2>(varrho, zero.data(), pn_pp_o.data());
			DR::Gen<2>(varrho, zero.data(), pn_lap_o.data());
		}

		template <>
		void init_x<2>() {
			varrho = 1./(1.*dp);
			Vec zero = Vec::Zero();
			DR::Gen<1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Gen<0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Gen<2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Gen<1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Gen<0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			pn_lap_o = pn_pp_o.block<1, PN::value>(0, 0) + pn_pp_o.block<1, PN::value>(2, 0);
		}

		template <>
		void init_x<3>() {
			Vec zero = Vec::Zero();
			DR::Gen<1, 0, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(0, 0).data());
			DR::Gen<0, 1, 0>(varrho, zero.data(), pn_p_o.block<1, PN::value>(1, 0).data());
			DR::Gen<0, 0, 1>(varrho, zero.data(), pn_p_o.block<1, PN::value>(2, 0).data());
			DR::Gen<2, 0, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(0, 0).data());
			DR::Gen<1, 1, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(1, 0).data());
			DR::Gen<1, 0, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(2, 0).data());
			DR::Gen<0, 2, 0>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(3, 0).data());
			DR::Gen<0, 1, 1>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(4, 0).data());
			DR::Gen<0, 0, 2>(varrho, zero.data(), pn_pp_o.block<1, PN::value>(5, 0).data());
			pn_lap_o = pn_pp_o.block<1, PN::value>(0, 0) + pn_pp_o.block<1, PN::value>(3, 0) + pn_pp_o.block<1, PN::value>(5, 0);
		}

		void init_x() { 
			invNeu.clear();
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(MatPP());
				if (type[p] == BD1) invNeu[p] = MatPP::Zero();
			}
			init_x<>();
		}

	public:
		std::vector<MatPP> invMat;
		std::unordered_map<unsigned, MatPP> invNeu;

		R varrho;
		Eigen::Matrix<R,D,PN::value,Eigen::RowMajor>					pn_p_o;
		Eigen::Matrix<R,mMath::H<D,2>::value,PN::value,Eigen::RowMajor>	pn_pp_o;
		Eigen::Matrix<R,1,PN::value,Eigen::RowMajor>					pn_lap_o;
	};

}