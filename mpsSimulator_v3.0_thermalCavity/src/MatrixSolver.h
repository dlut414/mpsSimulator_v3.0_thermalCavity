/*
*/
#pragma once
#include <iostream>
#include "Header.h"
#include "SparseMatrix.h"

namespace SIM {

	template <typename R, unsigned D>
	class MatrixSolver {
		typedef SparseMatrix<R> sMat;
		typedef DenseVector<R> dVec;
	public:
		MatrixSolver(const unsigned& _n, const R& e)
			: n(_n), au(D*_n, D*_n), u(D*_n), rhs(D*_n), eps(e) {
			init();
		}
		~MatrixSolver() {}

		void biCgStab_DiagonalPreconditioner() {
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void biCg_v() {
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void ccBiCg_augment(const std::vector<enum pType>& type) {
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}

	public:
		int n;
		int maxIter;
		R eps;
		sMat au;
		dVec u, rhs;

	private:
		void init() {
			maxIter = 1000;
			for (unsigned i = 0; i < n; i++) {
				x[i] = b[i] = 0.;
			}
			for (unsigned i = 0; i < D*n; i++) {
				u[i] = rhs[i] = 0.;
			}
		}
		void fina() {}
	};

}