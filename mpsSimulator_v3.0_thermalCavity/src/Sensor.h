/*
*/
#pragma once
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "Particle.h"

namespace SIM {

	template <typename R, unsigned D>
	class Sensor {
		typedef Eigen::Matrix<R,D,1> Vec;
	public:
		Sensor() {}
		~Sensor() {}

		void operator >> (const std::string& str) const {
			std::ofstream file("./out/s" + str + ".out", std::ofstream::out);
			for (auto s = 0; s < pos.size(); s++) {
				for (auto d = 0; d < D; d++) {
					file << std::setprecision(6) << std::scientific << pos[s][d] << " ";
				}
				for (auto d = 0; d < D; d++) {
					file << std::setprecision(6) << std::scientific << vect[s][d] << " ";
				}
				file << std::endl;
			}
			file.close();
			std::cout << " Writing Sensor. done " << std::endl;
		}
		void operator << (const std::string& str) {
			std::ifstream file(str);
			if (!file.is_open()) std::cout << " No file Sensor. file ! " << std::endl;
			pos.clear();	vect.clear();	scal.clear();
			while (file.good()) {
				Vec p;
				for (auto d = 0; d < D; d++) {
					file >> p[d];
				}
				pos.push_back(p);
				vect.push_back(Vec::Zero());
				scal.push_back(0.);
			}
			file.close();
			std::cout << " Reading Sensor. done " << std::endl;
		}

		void profile(const R& time, const std::string& str) const {
			std::ofstream file("./out/s" + str + ".out", std::ofstream::app);
			for (unsigned s = 0; s < pos.size(); s++) {
				file << std::setprecision(6) << std::scientific << time << " "
					<< std::setprecision(6) << std::scientific << scal[s]
					<< std::endl;
			}
			file.close();
			std::cout << " Writing profile. done " << std::endl;
		}

		template <typename R, unsigned D, typename Der>
		void writeVect(const Particle<R,D,Der>* const part) {
			const auto& pt = part->derived();
			for (auto s = 0; s < pos.size(); s++) {
				vect[s] = pt.func(pt.vel1, pos[s]);
			}
		}
		template <typename R, unsigned D, typename Der>
		void writeScal(const Particle<R, D, Der>* const part) {
			const auto& pt = part->derived();
			for (unsigned s = 0; s < pos.size(); s++) {
				scal[s] = pt.func(pt.pres, pos[s]);
			}
		}

	public:
		std::vector<Vec> pos;
		std::vector<Vec> vect;
		std::vector<R> scal;
	};

}