/*
 * LICENCE
 * copyright 2014 ~ ****
 * Some rights reserved.
 * Author: HUFANGYUAN
 * Released under CC BY-NC
*/
#pragma once
#include <cstdio>
#include "FractionalStep_x.h"
#include "FractionalStep_x_cst.h"
#include "Header.h"
#include "Controller.h"
#include "DrawSim.h"
#include "Bitmap.h"

struct Parameters {
	typedef double DataType;
	enum { Dim = 2, };
	enum { Order = 2, };
};

typedef SIM::FractionalStep_x_cst<Parameters::DataType, Parameters::Dim, Parameters::Order> Sim;
typedef REN::DrawSim<Parameters::DataType, Parameters::Dim> Ren;

Ren* renObj;
Sim* simObj;
REN::Controller stateObj;

namespace REN {

	void InitGL(int argc, char** argv);
	void MainLoop();
	void Final();

	void initOBJ();

	void Fps();
	void myMouse(int, int, int, int);
	void myMotion(int, int);
	void myMouseWheel(int, int, int, int);
	void myReshape(int, int);
	void myKeyboard(unsigned char, int, int);
	void myDisplay();

	void InitGL(int argc, char** argv) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
		glutInitWindowPosition(0, 0);
		glutInitWindowSize(stateObj.u_width, stateObj.u_height);
		glutCreateWindow("RTRenderer");
		//glLightModeli           (GL_LIGHT_MODEL_TWO_SIDE, 1);
		glutMouseFunc(myMouse);
		glutMotionFunc(myMotion);
		//glutMouseWheelFunc      (myMouseWheel);
		glutReshapeFunc(myReshape);
		glutKeyboardFunc(myKeyboard);
		glutDisplayFunc(myDisplay);

		glEnable(GL_TEXTURE_1D);
		glEnable(GL_TEXTURE_2D);
		glEnable(GL_TEXTURE_3D);
		glEnable(GL_CULL_FACE);
		//glDisable               (GL_CULL_FACE);
		glFrontFace(GL_CCW);
		glEnable(GL_POINT_SPRITE_ARB);
		glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);
		glEnable(GL_ALPHA_TEST);
		glAlphaFunc(GL_GREATER, 0.f);
		//glEnable                (GL_BLEND);
		//glBlendFunc             (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glClearColor(1.f, 1.f, 1.f, 0.f);

		glewInit();
	}
	
	void initOBJ() {
		simObj = new Sim;
		*simObj << "para.txt";
		simObj->init();
		renObj = new Ren(&stateObj);
		renObj->init();
	}

	void Fps() {}

	void MainLoop() { glutMainLoop(); }

	void Final() {}

	void callBack() {
		if (stateObj.b_save) {
			simObj->saveData();
			stateObj.b_save = false;
		}
		if (stateObj.b_sens) {
			simObj->sensorOut();
			stateObj.b_sens = false;
		}
		if (stateObj.b_bmp) {
			///write to bmp
			static Bitmap bm;
			static int i = 0;
			char name[256];
			sprintf_s(name, "./out/bm%04d.bmp", i++);
			bm.SaveAsBMP(name);
			stateObj.b_bmp = false;
		}
		if (!stateObj.b_stop) {
			double t = double(simObj->stepGL());
			static int count = int(t);
			if (t - 1.* count >= 0.) {
				stateObj.b_bmp = true;
				stateObj.b_sens = true;
				count++;
			}
			stateObj.b_dirty = true;
		}
	}

	void myMouse(int button, int s, int x, int y) {
		stateObj.clickMouse(button, s, x, y);
	}

	void myMotion(int x, int y) {
		stateObj.moveMouse(x, y);

		glutPostRedisplay();
	}

	void myMouseWheel(int button, int dir, int x, int y) {
		stateObj.rollMouse(button, dir, x, y);
	}

	void myReshape(int width, int height) {
		glViewport(0, 0, width, width);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluPerspective(-90.0f, float(stateObj.u_width) / float(stateObj.u_height), 1.0f, 100.0f);

		stateObj.reshapeWindow();
	}

	void myKeyboard(unsigned char key, int a, int b) {
		glutPostRedisplay();
		stateObj.pressKey(key, a, b);
		//if (stateObj.b_init) initOBJ();
		//stateObj.b_init = false;
		callBack();
	}

	void myDisplay() {
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), stateObj.m_pan)
			* glm::toMat4(stateObj.m_rotation)
			* glm::scale(glm::mat4(1.0f), stateObj.m_scale);

		stateObj.m_viewModelMat = stateObj.m_camera.GetViewMatrix() * modelMatrix;
		stateObj.m_projectionMat = stateObj.m_camera.GetProjectionMatrix();
		stateObj.m_projectionMatInv = glm::inverse(stateObj.m_projectionMat);
		stateObj.m_mvp = stateObj.m_projectionMat * stateObj.m_viewModelMat;
		stateObj.m_mvpInv = glm::inverse(stateObj.m_mvp);

		renObj->draw(simObj->part->type, simObj->part->pos, simObj->part->vort, simObj->part->pres);

		glutSwapBuffers();
		glutReportErrors();
		
		callBack();

		Fps();

		if (stateObj.b_dirty) {
			glutPostRedisplay();
			stateObj.b_dirty = false;
		}
		if (stateObj.b_leave) {
			glutLeaveMainLoop();
		}
	}

}

