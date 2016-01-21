/*
*/
// mpsSimulator_v2.0.cpp : Defines the entry point for the console application.
//

#include <tchar.h>
#include <omp.h>
#include <Eigen/Core>
#include "Header.h"
#include "renderer/Renderer.h"

int _tmain(int argc, _TCHAR* argv[]) {

	CreateDirectoryA(std::string(".\\out").c_str(), NULL);

	omp_set_num_threads(NTHREAD);
	Eigen::setNbThreads(NTHREAD);
	std::cout << "	Working Cores:	" << Eigen::nbThreads() << std::endl;

	REN::InitGL(argc, (char**)argv);
	REN::initOBJ();
	REN::MainLoop();
	REN::Final();

	return 0;
}

