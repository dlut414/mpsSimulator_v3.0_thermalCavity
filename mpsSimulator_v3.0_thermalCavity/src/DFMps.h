/*
*/
#pragma once
#include "Simulator.h"

namespace SIM {

	template <typename real, enum Dim dim>
	class DFMps : public Simulator< real, dim, DFMps<real, dim> > {
		typedef Vec3<real> vec;
		typedef Eigen::Triplet<real> tpl;
		typedef Eigen::SparseMatrix<real> spMat;
		typedef Eigen::Matrix<real, Eigen::Dynamic, 1> arr;
	public:
		DFMps() {}
		~DFMps() {}

		void step() {
			convect1();
			makeMat();
			makeSource();
			solvMat();
			convect2();
		}

	public:
		spMat a;
		arr x, b;

	private:
		void convect1() {
			for (unsigned i = 0; i < part->vel2.size(); i++) {
				if (part->type[i] != FLUID) continue;
				part->vel2[i] = part->vel1[i] + para.dt * (para.g + para.niu * part->lap(part->vel1, i));
			}
		}

		void convect2() {
			for (unsigned i = 0; i < part->vel2.size(); i++) {
				if (part->type[i] != FLUID) continue;
				switch (dim) {
				case TWOD:
					part->vel2[i].x += x[i*dim];
					part->vel2[i].z += x[i*dim+1];
					break;
				case THREED:
					part->vel2[i].x += x[i*dim];
					part->vel2[i].y += x[i*dim + 1];
					part->vel2[i].z += x[i*dim + 2];
					break;
				}
				part->pos[i] += 0.5 * (part->vel1[i] + part->vel2[i]);
				part->vel1[i] = part->vel2[i];
			}
		}

		void makeMat() {
			unsigned n = part->vel1.size();
			a = spMat(n, dim* n);
			x = arr(dim*n);	b = arr(n);
			std::vector<tpl> coef;
			coef.clear();
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) {
					for (int offset = 0; offset < dim; offset++) {
						coef.push_back(tpl(p, dim*p + offset, 1.));
					}
					continue;
				}
				vec pp = vec(0., 0., 0.);
				const iVec3 c = part->cell->iCoord(part->pos[p]);
				for (int k = -1; k <= 1; k++) {
					for (int j = -1; j <= 1; j++) {
						for (int i = -1; i <= 1; i++) {
							const iVec3 ne = c + iVec3(i, j, k);
							const unsigned key = part->cell->hash(ne);
							for (unsigned m = 0; m < part->cell->linkList[key].size(); m++) {
								const unsigned q = part->cell->linkList[key][m];
								if (q == p) continue;
								const vec	dr = part->pos[q] - part->pos[p];
								const real dr1 = dr.mag();
								const real dr2 = dr1 * dr1;
								const real w = part->w1(dr1);
								const vec pq = dr * (w / dr2);
								switch (dim) {
								case TWOD: 
									pp.x -= pq.x;
									pp.z -= pq.z;
									coef.push_back(tpl(p, q, pq.x));
									coef.push_back(tpl(p, q + 1, pq.z));
									break;
								case THREED: 
									pp.x -= pq.x;
									pp.y -= pq.y;
									pp.z -= pq.z;
									coef.push_back(tpl(p, q, pq.x));
									coef.push_back(tpl(p, q + 1, pq.y));
									coef.push_back(tpl(p, q + 2, pq.z));
									break;
								}
							}
						}
					}
				}
				switch (dim) {
				case TWOD:
					coef.push_back(tpl(p, p, pp.x));
					coef.push_back(tpl(p, p + 1, pp.z));
					break;
				case THREED:
					coef.push_back(tpl(p, p, pp.x));
					coef.push_back(tpl(p, p + 1, pp.y));
					coef.push_back(tpl(p, p + 2, pp.z));
					break;
				}
			}
			a.setFromTriplets(coef.begin(), coef.end());
		}

		void makeSource() {
			for (unsigned p = 0; p < part->pos.size(); p++) {
				if (part->type[p] != FLUID) {
					b[p] = 0.;
					continue;
				}
				b[p] = 0. - part->div(part->vel2, p);
			}
		}

		void solvMat() {
			spMat at = a.transpose();
			spMat aat = (a * at).pruned(para.eps);
			Eigen::SimplicialLDLT<spMat> solver;
			solver.compute(aat);
			unsigned n = part->vel1.size();
			spMat i(n, n);
			i.setIdentity();
			spMat inv = solver.solve(i);
			x = at * inv * mSol->b;
		}
	};

}