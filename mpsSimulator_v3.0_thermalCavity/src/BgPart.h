/*
*/
#pragma once
#include <vector>
#include "header.h"

namespace SIM {

	template <class real, enum Dim dim>
	class BgPart {
		typedef Vec3<real> vec;
		typedef Mat3<real> mat;
		typedef Eigen::Matrix<real, dim, dim> matEi;
	public:
		BgPart(const vec& p0, const real& dp, const real& _kr) {
			pos.clear();
			const int kr = int(_kr);
			switch (dim) {
			case TWOD:
				for (int k = -kr; k <= kr; k++) {
					for (int i = -kr; i <= kr; i++) {
						if (i == 0 && k == 0) continue;
						pos.push_back(vec(i*dp + p0.x, p0.y, k*dp + p0.z));
						counter.push_back(0);
					}
				}
				break;
			case THREED:
				for (int k = -kr; k <= kr; k++) {
					for (int j = -kr; j <= kr; j++) {
						for (int i = -kr; i <= kr; i++) {
							if (i == 0 && j == 0 && k == 0) continue;
							pos.push_back(vec(i*dp + p0.x, j*dp + p0.y, k*dp + p0.z));
							counter.push_back(0);
						}
					}
				}
				break;
			}
		}
		~BgPart() {}

	public:
		std::vector<vec> pos;
		std::vector<unsigned> counter;
	};

}