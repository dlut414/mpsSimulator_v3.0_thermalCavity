/*
*/
#pragma once
#include <vector>
#include <Eigen/Dense>
#include "Header.h"
#define _USE_MATH_DEFINES
#include <cmath>

namespace SIM {

	template <typename R, unsigned D>
	class ConsValue {
		typedef Eigen::Matrix<R, D, 1> vec;
		typedef Eigen::Matrix<R, 2, 1> vec2;
		typedef Eigen::Matrix<R, 3, 1> vec3;
		typedef Eigen::Matrix<R, 5, 1> vec5;
		typedef Eigen::Matrix<R, 6, 1> vec6;
		typedef Eigen::Matrix<R, 7, 1> vec7;
		typedef Eigen::Matrix<R, 8, 1> vec8;
		typedef Eigen::Matrix<R, 9, 1> vec9;
	public:
		ConsValue() {}
		~ConsValue() {}

		inline const R ww(const R& r, const R& re) const {
			if (r >= re) {
				return 0.;
			}
			else {
				return pow((1 - r / re), 2);
			}
		}

		inline const R w1(const R& r) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				return pow((1 - r / r0), 2);
			}
		}
		inline const R w1(const R& r, const R& re) const {
			if (r >= re) {
				return 0.;
			}
			else {
				return pow((1 - r / re), 2);
			}
		}

		inline const R w2(const R& r) const {
			if (r >= r0) {
				return 0.;
			}
			else {
				return pow((1 - r / r0), 2);
			}
		}
		inline const R w3(const R& r) const {
			if (r >= r0) return 0.;
			else {
				return pow((1. - r / r0), 2);
			}
		}
		inline const R w4_bSpline(const R& r) const {
			if (r >= r0) return 0.;
			const R q = r / r0;
			if (q <= 0.5) return (1. - 6 * q*q + 6 * q*q*q);
			else return (2.*pow((1. - q), 3));
		}

	public:
		R dp;
		R k;
		R r0;
		R beta;	//beta: free surface
		R n0;
		R lambda;
		unsigned pn0;
		R eps;
		R eps_mat;
		
	public:
		void init(const R& _k, const R& _beta) {
			k = _k; beta = _beta; r0 = k* dp;
			eps = std::numeric_limits<R>::epsilon();
			eps_mat = static_cast<R>(1.e-5);

			initConst();
		}

	private:
		void initConst() {
			std::vector<vec> v;
			if (D == 2) {
				v.clear();
				for (int i = 0; i < k * 2 + 1; i++) {
					for (int j = 0; j < k * 2 + 1; j++) {
						vec p;
						p << R(i), R(j);
						v.push_back(p);
					}
				}
				//cal n0
				//cal lambda
				n0 = 0.;
				lambda = 0.;
				pn0 = 0;
				for (unsigned i = 0; i < v.size(); i++) {
					R n = 0.;
					R lam = 0.;
					unsigned pn = 0;
					for (unsigned j = 0; j < v.size(); j++) {
						//if (j == i) continue;
						R w = w1((v[i] - v[j]).norm(), k);
						n += w;
						lam += (v[i] - v[j]).squaredNorm() * w;
						if (w > 0.) pn++;
					}
					lam = lam / n;
					n0 = n>n0 ? n : n0;
					lambda = lam > lambda ? lam : lambda;
					pn0 = pn > pn0 ? pn : pn0;
				}
			}
			else if (D == 3) {
				v.clear();
				for (int i = 0; i < k * 2 + 1; i++) {
					for (int j = 0; j < k * 2 + 1; j++) {
						for (int k = 0; k < k * 2 + 1; k++) {
							vec p;
							p << R(i), R(j), R(k);
							v.push_back(p);
						}
					}
				}
				//cal n0
				//cal lambda
				n0 = 0.;
				lambda = 0.;
				pn0 = 0;
				for (unsigned i = 0; i < v.size(); i++) {
					R n = 0.;
					R lam = 0.;
					unsigned pn = 0;
					for (unsigned j = 0; j < v.size(); j++) {
						//if (j == i) continue;
						R w = w1((v[i] - v[j]).norm(), k);
						n += w;
						lam += (v[i] - v[j]).squaredNorm() * w;
						if (w > 0.) pn++;
					}
					lam = lam / n;
					n0 = n>n0 ? n : n0;
					lambda = lam>lambda ? lam : lambda;
					pn0 = pn > pn0 ? pn : pn0;
				}
			}
			else {
				std::cout << " Wrong dimension number ! " << std::endl;
			}
			v.clear();
		}
	};

}