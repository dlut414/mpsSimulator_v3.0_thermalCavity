/*
*/
#pragma once

#define DEBUG 1
#define LEGACY 0
#define NTHREAD 4
#define OMP	1
#define BD_OPT 0

enum pType {
	FLUID = 0,
	BD1 = 1,
	BD2 = 2,
};
enum Dim {
	TWOD = 2,
	THREED = 3,
};
enum dType {
	SCALAR,
	VECTOR,
};

#define V_DIRICHLET		0x00000001
#define V_NEUMANN		0x00000002
#define P_DIRICHLET		0x00000004
#define P_NEUMANN		0x00000008
#define T_DIRICHLET0	0x00000010
#define T_DIRICHLET1	0x00000020
#define T_NEUMANN		0x00000040
#define T_DIRICHLET		(T_DIRICHLET0|T_DIRICHLET1)

#define ON(A, B)	(A|B)
#define OFF(A, B)	(A&(~B))
#define IS(A, B)	(A&B)