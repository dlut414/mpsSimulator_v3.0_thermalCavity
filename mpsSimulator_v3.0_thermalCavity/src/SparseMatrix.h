/*
* LICENCE
* copyright 2014 ~ ****
* Some rights reserved.
* Author: HUFANGYUAN
* Released under CC BY-NC
*/
#pragma once
#include <vector>

namespace SIM {

	template <typename R> class SparseMatrix;

	template <typename R>
	class DenseVector {
	public:
		DenseVector() : size(0) {}
		DenseVector(const int& n) :size(n) { value.resize(n); }
		template <typename U>
		DenseVector(const DenseVector<U>& dv) {
			size = dv.size;
			for (int i = 0; i < size; i++) {
				value.push_back(static_cast<R>(value[i]));
			}
		}
		~DenseVector() {}

		__forceinline R& operator[] (const int& i) {
			return value[i];
		}
		__forceinline const R& operator[] (const int& i) {
			return value[i];
		}
		__forceinline void copy(const DenseVector<R>& x) {
			for (int i = 0; i < size; i++) value[i] = x.value[i];
		}
		__forceinline void add(const DenseVector<R>& x) {
			for (int i = 0; i < size; i++) value[i] += x.value[i];
		}
		__forceinline void sub(const DenseVector<R>& x) {
			for (int i = 0; i < size; i++) value[i] -= x.value[i];
		}
		__forceinline void scale(const R& s) {
			for (int i = 0; i < size; i++) value[i] *= s;
		}
		__forceinline void axpy(const R& a, const DenseVector<R>& x) {
			for (int i = 0; i < size; i++) value[i] += a* x.value[i];
		}
		__forceinline const R dot(const DenseVector<R>& x) const {
			R sum = R(0);
			for (int i = 0; i < size; i++) sum += value[i] * x.value[i];
			return sum;
		}
		__forceinline void Ax(const SparseMatrix<R>& A, const DenseVector<R>& x) {
			for (int i = 0; i < size; i++) value[i] = R(0);
			for (int i = 0; i < A.size; i++) value[A.row[i]] += A.value[i] * x[A.col[i]];
		}

	public:
		int size;
		std::vector<R> value;
	};

	template <typename R>
	class SparseMatrix {
	public:
		SparseMatrix() : size(0) {}
		SparseMatrix(const int& nr, const int& nc) : size(0), nRow(nr), nCol(nc) {}
		~SparseMatrix() {}

		void clear() {
			row.clear();
			col.clear();
			value.clear();
			size = 0;
		}
		__forceinline void push(const int& r, const int& c, const R& v) {
			row.push_back(r);
			col.push_back(c);
			value.push_back(v);
			size++;
		}

	public:
		int size;
		int nRow;
		int nCol;
		std::vector<int> row;
		std::vector<int> col;
		std::vector<R> value;

	};

}