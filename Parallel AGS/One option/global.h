#pragma once
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <stdio.h>
#include <vector>
#include <queue>
#include <list>

#include <omp.h>

struct block
{
	double left;
	double right;
};

struct Pointer
{
	double x, z;
	int steps;
};
class global
{
private:

	long double time;
	double left = 0; // левая граница
	double right = 0; // правая граница
	double r = 0; //параметр для оценки константы Липшеца
	int steps=0;

	std::list<double> Values_Z;
	std::list<double> Arguments_X;

public:

	global();
	Pointer ParallelSolve(int threads, const double _left, const double _right, const double _r, const int _N_max, const double _Eps);
	Pointer solve(const double _left, const double _right, const double _r, const int _N_max, const double _Eps);
	inline double R(const double _m_small, const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double func(const double _x);
	double GetTime();
	int GetSteps() { return steps; };
};

#endif
