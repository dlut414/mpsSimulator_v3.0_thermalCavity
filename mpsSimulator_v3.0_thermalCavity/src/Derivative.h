/*
* LICENSE
*/
///Derivative.h
#pragma once
#include <vector>
#include "mMath.h"
#include "Polynomial.h"

namespace mMath {

	template <typename R, unsigned D, unsigned P>
	class Derivative_A_ {
	public:
		template <unsigned N>	static __forceinline R Power(const R& x)	{ return x*Power<N-1>(x); }
		template <>				static __forceinline R Power<0>(const R& x) { return 1; }
		template <unsigned N, unsigned M>	static __forceinline R A() { return (Factorial<N>::value) / (Factorial<N-M>::value); }
	};

	template <typename R, unsigned D, unsigned P>	class Derivative_A : public Derivative_A_<R,D,P> {};

	template <typename R, unsigned P>				class Derivative_A<R,1,P> : public Derivative_A_<R,1,P>{
	private:
		template <unsigned Der, unsigned P_ = 1, bool Zero = (Der>P_), bool Over = (P_==P)>
		struct Gen_ {
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = A<P_,Der>()*Power<P_-Der>(x); Gen_<Der,P_+1,(Der>(P_+1))>::Gen(x, out); }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = A<P_,Der>()*Power<P_-Der>(s*x)*Power<Der>(s); Gen_<Der,P_+1,(Der>(P_+1))>::Gen(s, x, out); }
		};
		template <unsigned Der, unsigned P_, bool Over>
		struct Gen_<Der,P_,true,Over> {
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = 0; Gen_<Der,P_+1,(Der>(P_+1))>::Gen(x, out); }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = 0; Gen_<Der,P_+1, (Der>(P_+1))>::Gen(s, x, out); }
		};
		template <unsigned Der, unsigned P_>
		struct Gen_<Der,P_,true,true> {
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = 0; }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = 0; }
		};
		template <unsigned Der, unsigned P_>
		struct Gen_<Der,P_,false,true> {
			static __forceinline void Gen(const R& x, R* const out)				{ out[P_-1] = A<P_,Der>()*Power<P_-Der>(x); }
			static __forceinline void Gen(const R& s, const R& x, R* const out) { out[P_-1] = A<P_,Der>()*Power<P_-Der>(s*x)*Power<Der>(s); }
		};
	public:
		template <unsigned Der> static __forceinline void Gen(const R& x, R* const out)				{ Gen_<Der>::Gen(x, out); }
		template <unsigned Der> static __forceinline void Gen(const R& s, const R& x, R* const out)	{ Gen_<Der>::Gen(s, x, out); }
	};

	template <typename R, unsigned P>				class Derivative_A<R,2,P> : public Derivative_A_<R,2,P> {
	private:
		template <unsigned DerX, unsigned DerY, unsigned P_ = 1, unsigned PX = P_, unsigned I = 0, bool Zero = (DerX>PX||DerY>(P_-PX))>
		struct Gen__									{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = A<PX,DerX>()*A<P_-PX,DerY>()*Power<PX-DerX>(x[0])*Power<P_-PX-DerY>(x[1]); Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = A<PX,DerX>()*A<P_-PX,DerY>()*Power<DerX+DerY>(s)*Power<PX-DerX>(s*x[0])*Power<P_-PX-DerY>(s*x[1]); Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned PX, unsigned I>
		struct Gen__ <DerX,DerY,P_,PX,I,true>			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = 0; Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = 0; Gen__<DerX,DerY,P_,PX-1,I+1>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen__<DerX,DerY,P_,0,I,true> 			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = 0; }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = 0; }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen__<DerX,DerY,P_,0,I,false> 			{
			static __forceinline void Gen(const R* const x, R* const out)				{ out[I] = A<0,DerX>()*A<P_-0,DerY>()*Power<0-DerX>(x[0])*Power<P_-0-DerY>(x[1]); }
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ out[I] = A<0,DerX>()*A<P_-0,DerY>()*Power<DerX+DerY>(s)*1*Power<P_-0-DerY>(s*x[1]); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_ = 1, unsigned I = 0, bool Over = (P_==P)>
		struct Gen_							{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<DerX,DerY,P_,P_,I>::Gen(x, out); Gen_<DerX,DerY,P_+1,I+H<2,P_>::value>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<DerX,DerY,P_,P_,I>::Gen(s, x, out); Gen_<DerX,DerY,P_+1,I+H<2,P_>::value>::Gen(s, x, out); }
		};
		template <unsigned DerX, unsigned DerY, unsigned P_, unsigned I>
		struct Gen_<DerX,DerY,P_,I,true>	{ 
			static __forceinline void Gen(const R* const x, R* const out)				{ Gen__<DerX,DerY,P_,P_,I>::Gen(x, out); } 
			static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen__<DerX,DerY,P_,P_,I>::Gen(s, x, out); }
		};
	public:
		template <unsigned DerX, unsigned DerY>	static __forceinline void Gen(const R* const x, R* const out)				{ Gen_<DerX,DerY>::Gen(x, out); }
		template <unsigned DerX, unsigned DerY>	static __forceinline void Gen(const R& s, const R* const x, R* const out)	{ Gen_<DerX,DerY>::Gen(s, x, out); }
	};
	


	template <typename R, unsigned D, unsigned P>
	class Derivative_B_ {
	public:
		enum { value = H<D,P>::value + Derivative_B_<R,D,P-1>::value, };
	};
	template <typename R, unsigned D>
	class Derivative_B_ <R,D,0> {
	public:
		enum { value = 1, };
	};

}