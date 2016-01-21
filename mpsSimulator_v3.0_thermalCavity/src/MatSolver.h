/*
*/
#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>
#include "header.h"

#define AUGMENT (1)
#define AG AUGMENT

namespace SIM {

	template <typename R, unsigned D>
	class MatSolver {
		typedef Eigen::Triplet<R> Tpl;
		typedef Eigen::Matrix<R,Eigen::Dynamic,1> dVec;
		typedef Eigen::SparseMatrix<R, Eigen::RowMajor> sMat;
#if AUGMENT
		typedef Eigen::IncompleteLUT<R> preconditioner;
#else
		typedef Eigen::DiagonalPreconditioner<R> preconditioner;
#endif
	public:
		MatSolver(const unsigned& _n, const R& e) 
			: n(_n), a(_n+AG,_n+AG), x(_n+AG), b(_n+AG), au(D*_n,D*_n), u(D*_n), rhs(D*_n), eps(e)
		{ 
			init();
		}
		~MatSolver() {}

		void biCg() {
			solverBiCg.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCg.solve(b);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void biCg_v() {
			solverBiCg.compute(au);
			//u = solver2.solveWithGuess(rhs, 0.5*u);
			u = solverBiCg.solve(rhs);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void qr() {
			a.makeCompressed();
			solverQR.compute(a);
			x = solverQR.solve(b);
			std::cout << " rank ----------------> " << solverQR.rank() << std::endl;
			if (solverQR.info() != Eigen::Success) {
				std::cout << " info ----------------> " << solverQR.lastErrorMessage() << std::endl;
			}
		}
		void lnBiCg() {
			sMat at(n, n);
			at = a.transpose();
			a = (a*at);
			solverBiCg.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCg.solve(b);
			x = (at*x);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void ccBiCg_augment(const std::vector<enum pType>& type) {
			sMat d(n + AG, n + AG);
			std::vector<Tpl> coef;
			for (unsigned p = 0; p<n; p++) {
				if (type[p] == BD2) continue;
				coef.push_back(Tpl(n, p, 1.));
				coef.push_back(Tpl(p, n, 1.));
			}
			d.setFromTriplets(coef.begin(), coef.end());
			a = a + d;
			b[n] = 0.;
			solverBiCg.compute(a);
			//x = solver1.solveWithGuess(b, 0.5*x);
			x = solverBiCg.solve(b);
			std::cout << " iterations ----------> " << solverBiCg.iterations() << std::endl;
			std::cout << " error ---------------> " << solverBiCg.error() << std::endl;
		}
		void ccCgs_horibata() {
			sMat at(n, n);
			at = a.transpose();
			dVec e(n), zero(n);
			solverBiCg.compute(at);
			e = solverBiCg.solve(0.*zero);
			const auto e2 = e.dot(e);
			if (e2 > eps) {
				b = b - (b.dot(e) / e2)*e;
			}
			const dVec r0 = b - a*x;
			auto r1 = r0;
			auto r2 = r1;
			auto p = r0;
			auto e1 = b;
			auto res = r1.dot(r1);
			auto h = b;
			int k;
			for (k = 0; k < maxIter; k++) {
				res = r1.dot(r1);
				if (res < eps) break;
				const auto alpha = res / r0.dot(a*p);
				h = e1 - alpha*(a*p);
				x += alpha*(e1 + h);
				r2 = r1 - alpha*(e1 + h);
				const auto beta = r0.dot(r2) / r0.dot(r1);
				e1 = r2 + beta* h;
				p = e1 + beta*(h + beta*p);
				r1 = r2;
			}
			std::cout << " iterations ----------> " << k << std::endl;
			std::cout << " error ---------------> " << res << std::endl;
		}
		void lsqr() {
		}

	public:
		unsigned n;
		int maxIter;
		R eps;
		sMat a, au;
		dVec x, b, u, rhs;
		Eigen::BiCGSTAB< sMat, preconditioner > solverBiCg;
		Eigen::SparseQR< sMat, Eigen::NaturalOrdering<int> > solverQR;

	private:
		void init() {
			maxIter = 1000;
			for (unsigned i = 0; i < n; i++) {
				x[i] = b[i] = 0.;
			}
			for (unsigned i = 0; i < D*n; i++) {
				u[i] = rhs[i] = 0.;
			}
#if AUGMENT
			solverBiCg.preconditioner().setDroptol(eps);
			solverBiCg.preconditioner().setFillfactor(1);
#endif
			solverBiCg.setMaxIterations(maxIter);
			solverBiCg.setTolerance(eps);
			solverQR.setPivotThreshold(1./n);
		}
		void fina() {}
	};

}