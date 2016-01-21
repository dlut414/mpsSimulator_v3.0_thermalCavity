/*
*/
#pragma once
#include "header.h"
#include "Particle.h"

#define UPWIND_VEL 0

namespace SIM {

	template <typename real, enum Dim dim>
	class Particle_2d_j : public Particle< real, dim, Particle_2d_j<real, dim> > {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, 1, 3> vec13;
		typedef Eigen::Matrix<real, 9, 1> vecp;
		typedef Eigen::Matrix<real, 9, 3> matp3;
		typedef Eigen::Matrix<real, 3, 3> mat33;
		typedef Eigen::Matrix<real, 9, 9> matpp;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		Particle_2d_j() : Particle() {}
		~Particle_2d_j() {}

		inline const vecp poly(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				s.x,
				s.z,
				s.x*s.x,
				s.x*s.z,
				s.z*s.z,
				s.x*s.x*s.x,
				s.x*s.x*s.z,
				s.x*s.z*s.z,
				s.z*s.z*s.z;
			return ret;
		}
		inline const vecp poly_px(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				1.,
				0.,
				2.*s.x,
				s.z,
				0.,
				3.*s.x*s.x,
				2.*s.z*s.x,
				s.z*s.z,
				0.;
			return ret / varrho[0];
		}
		inline const vecp poly_pz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				1.,
				0.,
				s.x,
				2.*s.z,
				0.,
				s.x*s.x,
				2.*s.x*s.z,
				3.*s.z*s.z;
			return ret / varrho[0];
		}
		inline const vecp poly_lap(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				2.,
				0.,
				2.,
				6.*s.x,
				2.*s.z,
				2.*s.x,
				6.*s.z;
			return ret / (varrho[0] * varrho[0]);
		}
		inline const vecp poly_pxx(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				2.,
				0.,
				0.,
				6.*s.x,
				2.*s.z,
				0.,
				0.;
			return ret / (varrho[0] * varrho[0]);
		}
		inline const vecp poly_pzz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				2.,
				0.,
				0.,
				2.*s.x,
				6.*s.z;
			return ret / (varrho[0] * varrho[0]);
		}
		inline const vecp poly_pxz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				1.,
				0.,
				0.,
				2.*s.x,
				2.*s.z,
				0.;
			return ret / (varrho[0] * varrho[0]);
		}
		inline const vecp poly_pxxx(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				0.,
				6.,
				0.,
				0.,
				0.;
			return ret / (varrho[0] * varrho[0] * varrho[0]);
		}
		inline const vecp poly_pxxz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				2.,
				0.,
				0.;
			return ret / (varrho[0] * varrho[0] * varrho[0]);
		}
		inline const vecp poly_pxzz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				2.,
				0.;
			return ret / (varrho[0] * varrho[0] * varrho[0]);
		}
		inline const vecp poly_pzzz(const vec& v) const {
			vecp ret;
			vec s = v / varrho[0];
			ret <<
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				0.,
				6.;
			return ret / (varrho[0] * varrho[0] * varrho[0]);
		}

		const real func(const std::vector<real>& phi, const unsigned& p) const {
			return phi[p];
		}

		const vec func(const std::vector<vec>& u, const unsigned& p) const {
			return u[p];
		}

		const vec grad(const std::vector<real>& phi, const unsigned& p) const {
			vecp  vv = vecp::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly(dr);
							vv += w * (phi[q] - phi[p]) * npq;
						}
					}
				}
			}
			vecp a = invMat[p] * vv;
			vecp px = poly_px_0;
			vecp pz = poly_pz_0;
			return vec(px.dot(a), 0., pz.dot(a));
		}

		const mat grad(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			auto a = invMat[p] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			return mat(vec(px.dot(a.block<9, 1>(0, 0)), 0., px.dot(a.block<5, 1>(0, 9))),
						vec(0.),
						vec(pz.dot(a.block<9, 1>(0, 0)), 0., pz.dot(a.block<5, 1>(0, 9))));
		}

		const real div(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly(dr);
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matp3 a = invMat[p] * vv;
			vecp px = poly_px_0;
			vecp pz = poly_pz_0;
			return a.block<9, 1>(0, 0).dot(px) + 0. + a.block<9, 1>(0, 2).dot(pz);
		}

		const real lap(const std::vector<real>& phi, const unsigned& p) const {
			vecp vv = vecp::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly(dr);
							vv += w * (phi[q] - phi[p]) * npq;
						}
					}
				}
			}
			vecp a = invMat[p] * vv;
			vecp lap = poly_lap_0;
			return lap.dot(a);
		}

		const vec lap(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const iVec3 c = cell->iCoord(pos[p]);
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec3 ne = c + iVec3(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const vec	dr = pos[q] - pos[p];
							const real	dr1 = dr.mag();
							if (dr1 > r0) continue;
							const real  w = w3(dr1);
							const vecp	npq = poly(dr);
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matp3 a = invMat[p] * vv;
			vecp lap = poly_lap_0;
			return vec(a.block<9, 1>(0, 0).dot(lap), 0., a.block<9, 1>(0, 2).dot(lap));
		}

		const vec rot(const std::vector<vec>& u, const unsigned& p) const {
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			auto a = invMat[p] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			return vec(0., px.dot(a.block<9, 1>(0, 2)) - pz.dot(a.block<9, 1>(0, 0)), 0.);
		}

		const real func(const std::vector<real>& phi, const vec& p) const {
			return real(0.);
		}

		const vec func(const std::vector<vec>& u, const vec& p) const {
			auto rid = 0;
			auto flag = 0;
			auto rr = std::numeric_limits<real>::max();
			auto c = cell->iCoord(p);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
							const auto dr = pos[q] - p;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							flag = 1;
							if (dr1 < rr) {
								rr = dr1;
								rid = q;
							}
						}
					}
				}
			}
			if (!flag) return vec(0.);

			matp3 vv = matp3::Zero();
			c = cell->iCoord(pos[rid]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(rid, q)) continue;
#endif
							const auto dr = pos[q] - pos[rid];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[rid].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[rid].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[rid].z) * npq;
						}
					}
				}
			}
			auto a = invMat[rid] * vv;
			auto px = poly_px_0;
			auto pz = poly_pz_0;
			auto gd = mat(vec(px.dot(a.block<9, 1>(0, 0)), 0., px.dot(a.block<9, 1>(0, 2))),
							vec(0.),
							vec(pz.dot(a.block<9, 1>(0, 0)), 0., pz.dot(a.block<9, 1>(0, 2))));
			return u[rid] + (p - pos[rid])*gd;
		}

		const vec func_nobd2(const std::vector<vec>& u, const vec& p) const {
			return  vec(0.);
		}

		const vec func_mafl(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = vec(0.);
				auto ww = real(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) ww = 1.;
				ret[fp] = ret[fp] / ww;
			}
			return ret[3] - (dmove.mag() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
		}

		const vec func_mafl_mmt(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto re = 1.5* dp;
			const auto dx = 1.5* dp;
			const auto p_i = pos[p];
			const auto dmove = p_new - p_i;
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dmove.norm();
#endif
			vec pLocal[6];

			for (int i = -3; i <= 2; i++) {
				pLocal[i + 3] = p_i - i*dx*up;
			}

			vec ret[6];
			ret[3] = u[p];
			for (int fp = 1; fp <= 4; fp++) {
				if (fp == 3) continue;
				ret[fp] = vec(0.);
				auto ww = real(0.);
				auto c = cell->iCoord(pLocal[fp]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
								const unsigned q = cell->linkList[key][m];
								if (type[q] == BD2) continue;
#if BD_OPT
								if (bdOpt(q)) continue;
#endif
								const auto dr1 = (pos[q] - pLocal[fp]).mag();
								const auto dr1_m1 = (pos[q] - pLocal[fp - 1]).mag();
								const auto dr1_p1 = (pos[q] - pLocal[fp + 1]).mag();
								if (dr1 > re) continue;
								if (dr1 > dr1_m1 || dr1 > dr1_p1) continue;
								const auto w = w1(dr1);
								ww += w;
								ret[fp] += w * u[q];
							}
						}
					}
				}
				if (abs(ww) < eps) continue;
				ret[fp] = ret[fp] / ww;
			}
			auto ret_mmt = ret[3] - (dmove.mag() / dx)* (0.125* ret[1] - 0.875* ret[2] + 0.375* ret[3] + 0.375* ret[4]);
			auto ret_min = ret[1];
			auto ret_max = ret[1];
			for (int fp = 1; fp <= 4; fp++) {
				if (ret[fp].mag2() < ret_min.mag2()) ret_min = ret[fp];
				if (ret[fp].mag2() > ret_max.mag2()) ret_max = ret[fp];
			}
			if (ret_mmt.mag2() < ret_min.mag2()) ret_mmt = ret_min* ret_mmt.norm();
			if (ret_mmt.mag2() > ret_max.mag2()) ret_mmt = ret_max* ret_mmt.norm();
			return ret_mmt;
		}

		const vec func_mls_a(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 9>(0, 0) += w * npq[0] * npq;
							mm.block<1, 9>(1, 0) += w * npq[1] * npq;
							mm.block<1, 9>(2, 0) += w * npq[2] * npq;
							mm.block<1, 9>(3, 0) += w * npq[3] * npq;
							mm.block<1, 9>(4, 0) += w * npq[4] * npq;
							mm.block<1, 9>(5, 0) += w * npq[5] * npq;
							mm.block<1, 9>(6, 0) += w * npq[6] * npq;
							mm.block<1, 9>(7, 0) += w * npq[7] * npq;
							mm.block<1, 9>(8, 0) += w * npq[8] * npq;
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<5, 5>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					auto mm__ = mm.block<2, 2>(0, 0);
					if (abs(mm__.determinant()) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv.block<2, 2>(0, 0) = mm__.inverse();
				}
				else inv.block<5, 5>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<9, 1>(0, 0);
			const auto az = a.block<9, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			const auto pxxx = poly_pxxx_0;
			const auto pxxz = poly_pxxz_0;
			const auto pxzz = poly_pxzz_0;
			const auto pzzz = poly_pzzz_0;
			const auto dp = p_new - pos[p];
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax))
				+ (1./6.)* ( dp.x*dp.x*dp.x*pxxx.dot(ax) + 3.*dp.x*dp.x*dp.z*pxxz.dot(ax) + 3.*dp.x*dp.z*dp.z*pxzz.dot(ax) + dp.z*dp.z*dp.z*pzzz.dot(ax) );
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az))
				+ (1. / 6.)* (dp.x*dp.x*dp.x*pxxx.dot(az) + 3.*dp.x*dp.x*dp.z*pxxz.dot(az) + 3.*dp.x*dp.z*dp.z*pxzz.dot(az) + dp.z*dp.z*dp.z*pzzz.dot(az));
			return ret;
		}

		const vec func_mls_a_upwind(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {
			const auto dp = p_new - pos[p];
#if UPWIND_VEL
			const auto up = -(u[p].norm());
#else
			const auto up = dp.norm();
#endif
			matpp mm = matpp::Zero();
			matp3 vv = matp3::Zero();
			const auto c = cell->iCoord(pos[p]);
			for (auto k = -1; k <= 1; k++) {
				for (auto j = -1; j <= 1; j++) {
					for (auto i = -1; i <= 1; i++) {
						const auto ne = c + iVec3(i, j, k);
						const auto key = cell->hash(ne);
						for (auto m = 0; m < cell->linkList[key].size(); m++) {
							const auto q = cell->linkList[key][m];
#if BD_OPT
							if (bdOpt(p, q)) continue;
#endif
							const auto dr = pos[q] - pos[p];
							if (dr*up < 0.) continue;
							const auto dr1 = dr.mag();
							if (dr1 > r0) continue;
							const auto w = w3(dr1);
							const auto npq = poly(dr);
							mm.block<1, 9>(0, 0) += w * npq[0] * npq;
							mm.block<1, 9>(1, 0) += w * npq[1] * npq;
							mm.block<1, 9>(2, 0) += w * npq[2] * npq;
							mm.block<1, 9>(3, 0) += w * npq[3] * npq;
							mm.block<1, 9>(4, 0) += w * npq[4] * npq;
							mm.block<1, 9>(5, 0) += w * npq[5] * npq;
							mm.block<1, 9>(6, 0) += w * npq[6] * npq;
							mm.block<1, 9>(7, 0) += w * npq[7] * npq;
							mm.block<1, 9>(8, 0) += w * npq[8] * npq;
							vv.block<9, 1>(0, 0) += w * (u[q].x - u[p].x) * npq;
							vv.block<9, 1>(0, 1) += w * (u[q].y - u[p].y) * npq;
							vv.block<9, 1>(0, 2) += w * (u[q].z - u[p].z) * npq;
						}
					}
				}
			}
			matpp inv = matpp::Zero();
			if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
				std::cout << mm.determinant() << std::endl;
#endif
				auto mm_ = mm.block<5, 5>(0, 0);
				if (abs(mm_.determinant()) < eps_mat) {
					auto mm__ = mm.block<2, 2>(0, 0);
					if (abs(mm__.determinant()) < eps_mat) {
						inv = matpp::Zero();
					}
					else inv.block<2, 2>(0, 0) = mm__.inverse();
				}
				else inv.block<5, 5>(0, 0) = mm_.inverse();
			}
			else inv = mm.inverse();

			const auto a = inv * vv;
			const auto ax = a.block<9, 1>(0, 0);
			const auto az = a.block<9, 1>(0, 2);
			const auto px = poly_px_0;
			const auto pz = poly_pz_0;
			const auto pxx = poly_pxx_0;
			const auto pzz = poly_pzz_0;
			const auto pxz = poly_pxz_0;
			const auto pxxx = poly_pxxx_0;
			const auto pxxz = poly_pxxz_0;
			const auto pxzz = poly_pxzz_0;
			const auto pzzz = poly_pzzz_0;
			auto ret = u[p];
			ret.x += dp.x* (px.dot(ax)) + dp.z* (pz.dot(ax))
				+ 0.5* dp.x* dp.x* (pxx.dot(ax)) + dp.x* dp.z* (pxz.dot(ax)) + 0.5* dp.z* dp.z* (pzz.dot(ax))
				+ (1. / 6.)* (dp.x*dp.x*dp.x*pxxx.dot(ax) + 3.*dp.x*dp.x*dp.z*pxxz.dot(ax) + 3.*dp.x*dp.z*dp.z*pxzz.dot(ax) + dp.z*dp.z*dp.z*pzzz.dot(ax));
			ret.z += dp.x* (px.dot(az)) + dp.z* (pz.dot(az))
				+ 0.5* dp.x* dp.x* (pxx.dot(az)) + dp.x* dp.z* (pxz.dot(az)) + 0.5* dp.z* dp.z* (pzz.dot(az))
				+ (1. / 6.)* (dp.x*dp.x*dp.x*pxxx.dot(az) + 3.*dp.x*dp.x*dp.z*pxxz.dot(az) + 3.*dp.x*dp.z*dp.z*pxzz.dot(az) + dp.z*dp.z*dp.z*pzzz.dot(az));
			return ret;
		}

		const vec func_mls_b(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {}

		const vec func_mls_b_upwind(const std::vector<vec>& u, const unsigned& p, const vec& p_new) const {}

		void updateInvMat() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] == BD2) continue;
				matpp mm = matpp::Zero();
				const auto c = cell->iCoord(pos[p]);
				for (auto k = -1; k <= 1; k++) {
					for (auto j = -1; j <= 1; j++) {
						for (auto i = -1; i <= 1; i++) {
							const auto ne = c + iVec3(i, j, k);
							const auto key = cell->hash(ne);
							for (auto m = 0; m < cell->linkList[key].size(); m++) {
								const auto q = cell->linkList[key][m];
#if BD_OPT
								if (bdOpt(p, q)) continue;
#endif
								const auto dr = pos[q] - pos[p];
								const auto dr1 = dr.mag();
								if (dr1 > r0) continue;
								const auto w = w3(dr1);
								const auto npq = poly(dr);
								mm.block<1, 9>(0, 0) += w * npq[0] * npq;
								mm.block<1, 9>(1, 0) += w * npq[1] * npq;
								mm.block<1, 9>(2, 0) += w * npq[2] * npq;
								mm.block<1, 9>(3, 0) += w * npq[3] * npq;
								mm.block<1, 9>(4, 0) += w * npq[4] * npq;
								mm.block<1, 9>(5, 0) += w * npq[5] * npq;
								mm.block<1, 9>(6, 0) += w * npq[6] * npq;
								mm.block<1, 9>(7, 0) += w * npq[7] * npq;
								mm.block<1, 9>(8, 0) += w * npq[8] * npq;
							}
						}
					}
				}
				invMat[p] = matpp::Zero();
				if (abs(mm.determinant()) < eps_mat) {
#if DEBUG
					std::cout << mm.determinant() << std::endl;
#endif
					auto mm_ = mm.block<5, 5>(0, 0);
					if (abs(mm_.determinant()) < eps_mat) {
						auto mm__ = mm.block<2, 2>(0, 0);
						if (abs(mm__.determinant()) < eps_mat) {
							invMat[p] = matpp::Zero();
							continue;
						}
						invMat[p].block<2, 2>(0, 0) = mm__.inverse();
						continue;
					}
					invMat[p].block<5, 5>(0, 0) = mm_.inverse();
					continue;
				}
				invMat[p] = mm.inverse();
			}
		}

		void updateDiver() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p<int(np); p++) {
				diver[p] = div(part->vel1, p);
			}
		}

		void init2d_x() {
			invMat.clear();
			for (int p = 0; p < int(np); p++) {
				invMat.push_back(matpp());
			}

			varrho.clear();
			varrho.push_back(1.*dp);
			poly_px_0 = poly_px(vec(0.));
			poly_pz_0 = poly_pz(vec(0.));
			poly_lap_0 = poly_lap(vec(0.));
			poly_pxx_0 = poly_pxx(vec(0.));
			poly_pzz_0 = poly_pzz(vec(0.));
			poly_pxz_0 = poly_pxz(vec(0.));
			poly_pxxx_0 = poly_pxxx(vec(0.));
			poly_pxxz_0 = poly_pxxz(vec(0.));
			poly_pxzz_0 = poly_pxzz(vec(0.));
			poly_pzzz_0 = poly_pzzz(vec(0.));
		}

	public:
		std::vector<matpp> invMat;

		std::vector<real> varrho;
		vecp poly_px_0;
		vecp poly_pz_0;
		vecp poly_lap_0;
		vecp poly_pxx_0;
		vecp poly_pzz_0;
		vecp poly_pxz_0;
		vecp poly_pxxx_0;
		vecp poly_pxxz_0;
		vecp poly_pxzz_0;
		vecp poly_pzzz_0;
	};

}