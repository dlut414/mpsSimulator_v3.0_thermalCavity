/*
*/
#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include "Header.h"
#include "LinkCell.h"
#include "ConsValue.h"

namespace SIM {
	
	template <typename R, unsigned D, typename Derived>
	class Particle : public ConsValue<R, D> {
		typedef Eigen::Matrix<int, D, 1>	iVec;
		typedef Eigen::Matrix<R, D, 1> vec;
		typedef Eigen::Matrix<R, 1, 3> vec13;
		typedef Eigen::Matrix<R, 5, 1> vec5;
		typedef Eigen::Matrix<R, 6, 1> vec6;
		typedef Eigen::Matrix<R, 7, 1> vec7;
		typedef Eigen::Matrix<R, 8, 1> vec8;
		typedef Eigen::Matrix<R, 9, 1> vec9;
		typedef Eigen::Matrix<R, 5, 5> mat55;
		typedef Eigen::Matrix<R, 6, 6> mat66;
		typedef Eigen::Matrix<R, 7, 7> mat77;
		typedef Eigen::Matrix<R, 8, 8> mat88;
		typedef Eigen::Matrix<R, 9, 9> mat99;
		typedef Eigen::Matrix<R, 5, 3> mat53;
		typedef Eigen::Matrix<R, 6, 3> mat63;
		typedef Eigen::Matrix<R, 7, 3> mat73;
		typedef Eigen::Matrix<R, 8, 3> mat83;
		typedef Eigen::Matrix<R, 9, 3> mat93;
		typedef Eigen::Matrix<R, D, D> mat;
	public:
		Particle() {}
		~Particle() {}

		Derived& derived() { return *static_cast<Derived*>(this); }
		const Derived& derived() const { return *static_cast<const Derived*>(this); }

		void clean() {
			type.clear(); pos.clear(); pos_m1.clear(); vel1.clear(); vel2.clear(); vel_m1.clear();
			temp.clear(); pnd.clear(); pres.clear(); pn.clear(); nbd.clear(); fs.clear();
			team.clear(); phi.clear(); vort.clear(); norm.clear(); bdc.clear();
		}
		void operator >> (const std::string str) const {
			std::ofstream file(str, std::ofstream::out);
			file << std::scientific << std::setprecision(6) << ct << std::endl;
			file << std::scientific << std::setprecision(6) << dp << std::endl;
			file << np << " " << bd1 << " " << bd2 << std::endl;
			for (unsigned p = 0; p < np; p++) {
				file << std::scientific << std::setprecision(6);
				file << type[p] << " ";
				for (int d = 0; d < D; d++) {
					file << pos[p][d] << " ";
				}
				for (int d = 0; d < D; d++) {
					file << vel1[p][d] << " ";
				}
				file << std::endl;
			}
			std::cout << " Writing Geo. done. " << std::endl;
			file.close();
		}
		void operator << (const std::string str) {
			int n;	 int t;		vec p;		vec	v;
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " File Geo. not found ! " << std::endl;
			file >> ct >> dp >> np >> bd1 >> bd2;
			n = np;
			while (n-- > 0) {
				file >> t;
				for (int d = 0; d < D; d++) file >> p[d];
				for (int d = 0; d < D; d++) file >> v[d];
				addPart(pType(t), p, v);
			}
			file.close();
			std::cout << " Reading Geo. done " << std::endl;
		}

		void addPart (const pType& t, const vec& p, const vec& v) {
			type.push_back(t);	pos.push_back(p); pos_m1.push_back(p);
			vel1.push_back(v);	vel2.push_back(v); vel_m1.push_back(v);
			temp.push_back(0.); pnd.push_back(0.);	pres.push_back(0.);
			pn.push_back(0);	nbd.push_back(0);
			fs.push_back(0);	team.push_back(0);	
			phi.push_back(0.); vort.push_back(0.);
			norm.push_back(vec::Zero()); bdc.push_back(0);
		}

		const R cPnd(const unsigned& p) const {
			R ret = 0.;
			const iVec c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					//if (q == p) continue;
					//if (type[q]==FLUID && (team[p]!=team[q])) continue;
					const R dr1 = (pos[q] - pos[p]).norm();
					ret += w1(dr1);
				}
			}
			return ret;
		}

		const R cPnd(const vec& p) const {
			R ret = 0.;
			const iVec c = cell->iCoord(p);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					const R dr1 = (pos[q] - p).norm();
					ret += w1(dr1);
				}
			}
			return ret;
		}

		const R cPn(const unsigned& p) const {
			R ret = 0;
			const iVec c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					if (q == p) continue;
					const R	dr1 = (pos[q] - pos[p]).norm();
					//if (dr1 < r0) ret += (isFs(q) ? 0.7 : 1.);
					if (dr1 < r0) ret += 1.;
				}
			}
			return ret;
		}

		const unsigned cNbd(const unsigned& p) const {
			unsigned ret = 0;
			const iVec c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					const R	dr1 = (pos[q] - pos[p]).norm();
					if (q == p || dr1 > r0) continue;
					if (isFs(q)) ret++;
				}
			}
			return ret;
		}

		void buildCell() {
			BBox<R> b = BBox<R>();
			for (unsigned p = 0; p < pos.size(); p++) {
				b += pos[p];
			}
			b.Expand(0.1);
			cell = new LinkCell<R,D>(b, r0);
			updateCell();
		}
		inline void updateCell() {
			cell->update(pos);
		}

		void updateTeam() {
			for (unsigned p = 0; p < np; p++) team[p] = 0;
			auto t = 1;
			for (unsigned p = 0; p < np; p++) {
				if (type[p] != FLUID) continue;
				dfs(p, t++);
			}
		}

		void dfs(const unsigned& p, const int& t) {
			if (type[p] == BD2) return;
			if (team[p] != 0) return;
			if (type[p] == BD1) { team[p] = t; return; }
			team[p] = t;
			const iVec c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					const R	dr1 = (pos[q] - pos[p]).norm();
					if (q == p || dr1 > 2.*dp) continue;
					dfs(q, t);
				}
			}
			return;
		}

		void b2b() {
			bbMap.clear();
			for (unsigned p = 0; p < pos.size(); p++) {
				if (type[p] != BD2) continue;
				auto tmpdr = std::numeric_limits<R>::max();
				unsigned tmpbb = 0;
				const auto c = cell->iCoord(pos[p]);
				for (auto i = 0; i < cell->blockSize::value; i++) {
					const auto key = cell->hash(c, i);
					for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
						const unsigned q = cell->linkList[key][m];
						if (q == p || type[q] != BD1) continue;
						const auto dr1 = (pos[q] - pos[p]).norm();
						if (dr1 < tmpdr) {
							tmpdr = dr1;
							tmpbb = q;
						}
					}
				}
				bbMap[p] = tmpbb;
			}
		}

		void b2norm() {
			bdnorm.clear();
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				const auto n = pos[p.second] - pos[p.first];
				const auto tmp = bdnorm.at(p.second);
				if(tmp.norm() < n.norm()) bdnorm[p.second] = n;
			}
			for (const auto& p : bbMap) {
				bdnorm[p.second] = bdnorm.at(p.second).normalized();
			}
		}

		void b2neumann() {
			p_neumann.clear();
			t_neumann.clear();
			for (unsigned p = 0; p < np; p++) {
				if (IS(bdc[p], P_NEUMANN)) p_neumann[p] = 0.;
				if (IS(bdc[p], T_NEUMANN)) t_neumann[p] = 0.;
			}
		}

		void b2dirichlet() {
			p_dirichlet.clear();
			t_dirichlet.clear();
			for (unsigned p = 0; p < np; p++) {
				if (IS(bdc[p], P_DIRICHLET)) p_dirichlet[p] = 0.;
				if (IS(bdc[p], T_DIRICHLET0)) t_dirichlet[p] = 0.;
				if (IS(bdc[p], T_DIRICHLET1)) t_dirichlet[p] = 1.;
			}
		}

		void makeBdc() {
			for (unsigned p = 0; p < np; p++) {
				bdc[p] = 0;
				if (type[p] == BD1) {
					bdc[p] = ON(bdc[p], P_NEUMANN);
					if (abs(pos[p][0] - 0.) < eps) bdc[p] = ON(bdc[p], T_DIRICHLET1);
					if (abs(pos[p][0] - 1.) < eps) bdc[p] = ON(bdc[p], T_DIRICHLET0);
					if (abs(pos[p][1] - 0.) < eps) bdc[p] = ON(bdc[p], T_NEUMANN);
					if (abs(pos[p][1] - 1.) < eps) bdc[p] = ON(bdc[p], T_NEUMANN);
				}
			}
		}

		inline const int isFs(const unsigned& p) const {
			return fs[p];
			//return _isFs(p);
		}

		inline const int isFs(const vec& p) const {
			return (cPnd(p) < beta * n0);
		}

		const int _isFs(const unsigned& p) {
			//return (pnd[p] < beta * n0);

			/*tamai surface detection*/
			if (type[p] == BD2) return 0;
			if (pnd[p] > (beta * n0)) return 0;
			mat mm = mat(0.);
			const iVec c = cell->iCoord(pos[p]);
			for (auto i = 0; i < cell->blockSize::value; i++) {
				const auto key = cell->hash(c, i);
				for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
					const unsigned q = cell->linkList[key][m];
					if (q == p) continue;
					const vec	dr = pos[q] - pos[p];
					const R	dr1 = dr.norm();
					if (dr1 > r0) continue;
					const R  w = w2(dr1);
					const vec	npq = dr / dr1;
					mm.x += w * npq.x * npq;
					mm.y += w * npq.y * npq;
					mm.z += w * npq.z * npq;
				}
			}
			mm = (D / n0) * mm;

			Eigen::SelfAdjointEigenSolver<mat> sol(mm);
			vec eig = sol.eigenvalues();
			mat eigvs = sol.eigenvectors();
			R eigmin = eig[0];
			vec eigv = eigvs.col(0);
			vec neigv = -eigv;

			if (eigmin <= 0.2) return 2;
			if (eigmin > 0.8) return 0;

			const R root2 = 1.415;
			int flag1 = 1, flag2 = 1;
			for (int k = -1; k <= 1; k++) {
				for (int j = -1; j <= 1; j++) {
					for (int i = -1; i <= 1; i++) {
						const iVec ne = c + iVec(i, j, k);
						const unsigned key = cell->hash(ne);
						for (unsigned m = 0; m < cell->linkList[key].size(); m++) {
							const unsigned q = cell->linkList[key][m];
							if (q == p) continue;
							const vec	dr = pos[q] - pos[p];
							const R	dr1 = dr.norm();
							if (dr1 < root2 * dp) {
								if ((dr / dr1) * eigv >(root2 / 2.)) flag1 = 0;
								if ((dr / dr1) * neigv >(root2 / 2.)) flag2 = 0;
							}
							else {
								if ((pos[p] + dp * eigv - pos[q]).mag() < dp) flag1 = 0;
								if ((pos[p] + dp * neigv - pos[q]).mag() < dp) flag2 = 0;
							}
						}
					}
				}
			}
			if (flag1) {
				norm[p] = eigv;
				return 1;
			}
			if (flag2) {
				norm[p] = neigv;
				return 1;
			}
			return 0;
		}

		inline const int bdOpt(const unsigned& p, const unsigned& q) const {
			//if (type[q] == BD2) return 1;
			//if (team[p] != team[q]) return 1;
			if (type[q] == BD1 && isFs(q)) return 1;
			//if (type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			//if (type[p] == FLUID && type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			return 0;
		}
		inline const int bdOpt(const unsigned& q) const {
			//if (type[q] == BD2) return 1;
			//if (team[p] != team[q]) return 1;
			if (type[q] == BD1 && isFs(q)) return 1;
			//if (type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			//if (type[p] == FLUID && type[q] == BD2 && isFs(bbMap.at(q))) return 1;
			return 0;
		}
		inline const int bdSlip(const unsigned& q) const {
			///slip
			if (type[q] != FLUID) return 1;
			return 0;
		}

		void bdNoSlip() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] != BD2) continue;
				const auto mid = bbMap[p];
				const auto mirror = 2.* pos[mid] - pos[p];
				const auto vmir = derived().func_nobd2(vel1, mirror);
				///axis symmetric
				//const auto norm = (pos[mid] - pos[p]).norm();
				//const auto vnor = (vmir* norm)* norm;
				//vel1[p] = vel2[p] = 2.* vnor - vmir;
				///point symmetric
				vel1[p] = vel2[p] = -vmir;
			}
		}
		void bdSetZero() {
#if OMP
#pragma omp parallel for
#endif
			for (int p = 0; p < int(np); p++) {
				if (type[p] != BD2) continue;
				vel1[p] = vel2[p] = vec(0.);
			}
		}

	public:
		R ct;
		unsigned np, bd1, bd2;
		std::vector<vec> pos;
		std::vector<vec> pos_m1;
		std::vector<vec> vel1;
		std::vector<vec> vel2;
		std::vector<vec> vel_m1;
		std::vector<vec> norm;

		std::vector<R>	temp;
		std::vector<R>	pnd;
		std::vector<R>	pres;
		std::vector<pType>	type;
		std::vector<int> bdc;
		std::vector<R>	pn;
		std::vector<unsigned> nbd;
		std::vector<int> fs;
		std::vector<int> team;
		std::vector<R> phi;
		std::vector<R> vort;
		std::unordered_map<unsigned, vec> bdnorm;
		std::unordered_map<unsigned, R> p_dirichlet;
		std::unordered_map<unsigned, R> t_dirichlet;
		std::unordered_map<unsigned, R> p_neumann;
		std::unordered_map<unsigned, R> t_neumann;
		std::unordered_map<unsigned, unsigned> bbMap;

		LinkCell<R,D>* cell;

	private:
	};

}