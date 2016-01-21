/*
* LICENSE
*/
///Polynomial.h
#pragma once
#include <vector>
#include "mMath.h"

namespace mMath {

	template <typename R, unsigned D, unsigned P>
	class Polynomial_A_ {
	public:
		enum { value = H<D, P>::value + Polynomial_A_<R,D,P-1>::value, };
		template <unsigned N>	static __forceinline R Power(const R& x) { return x*Power<N - 1>(x); }
		template <>				static __forceinline R Power<0>(const R& x) { return 1; }
	};
	template <typename R, unsigned D>	class Polynomial_A_<R,D,0> { public: enum { value = 0, }; };

	template <typename R, unsigned D, unsigned P>	class Polynomial_A : public Polynomial_A_<R,D,P>{};

	template <typename R, unsigned P>				class Polynomial_A<R,1,P> : public Polynomial_A_<R,1,P> {
	private:
		template <unsigned P_ = 1>	static __forceinline void Gen_		(const R& x, R* const out) { out[P_-1] = Power<P_>(x); Gen_<P_+1>(x, out); }
		template <>					static __forceinline void Gen_<P>	(const R& x, R* const out) { out[P-1] = Power<P>(x); }
		template <unsigned P_ = 1>	static __forceinline void Gen_		(const R& s, const R& x, R* const out) { out[P_-1] = Power<P_>(s*x); Gen_<P_+1>(s, x, out); }
		template <>					static __forceinline void Gen_<P>	(const R& s, const R& x, R* const out) { out[P-1] = Power<P>(s*x); }
	public:
		static __forceinline void Gen(const R& x, R* const out)				{ Gen_(x, out); }
		static __forceinline void Gen(const R& s, const R& x, R* const out) { Gen_(s, x, out); }
	};

	template <typename R, unsigned P>				class Polynomial_A<R,2,P> : public Polynomial_A_<R,2,P> {
	private:
		template <unsigned P_ = 1, unsigned PX = P_, unsigned I = 0>
		struct Gen__			{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = Power<PX>(x[0])*Power<P_-PX>(x[1]); Gen__<P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = Power<PX>(s*x[0])*Power<P_-PX>(s*x[1]); Gen__<P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned P_, unsigned I>						
		struct Gen__<P_,0,I>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = Power<P_>(x[1]); }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = Power<P_>(s*x[1]); }
		};
		template <unsigned P_ = 1, unsigned I = 0, bool Over = (P==P_)>
		struct Gen_				{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<P_,P_,I>::Gen(x, out); Gen_<P_+1,I+H<2,P_>::value>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<P_,P_,I>::Gen(s, x, out); Gen_<P_+1,I+H<2,P_>::value>::Gen(s, x, out); }
		};
		template <unsigned P_, unsigned I>
		struct Gen_<P_,I,true>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<P_,P_,I>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<P_,P_,I>::Gen(s, x, out); }
		};
	public:
		static __forceinline void Gen(const R* const x, R* const out)				{ Gen_<>::Gen(x, out); }
		static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen_<>::Gen(s, x, out); }
	};
	


	template <typename R, unsigned D, unsigned P>
	class Polynomial_B_ {
	public:
		enum { value = H<D, P>::value + Polynomial_B_<R, D, P-1>::value, };
	};
	template <typename R, unsigned D>
	class Polynomial_B_ < R, D, 0 > {
	public:
		enum { value = 1, };
	};

}