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

class global
{
private:
	double left = 0; // левая граница
	double right = 0; // правая граница
	double r = 0; // оценка константы Липшица

	std::list<block> z_block;
	std::list<block> x_block;
	std::list<block> M_block;
	std::list<block> R_block;

	std::list<double> Values_Z;
	std::list<double> Arguments_X;

	double solutionX = 0;
	double solutionZ = 0;
	int procs = 1;
	long double time;
	__int64 steps;

public:
	global();
	double solve(const double _left, const double _right, const double _r, const int _N_max, const double _Eps);
	inline double R(const double _m_small, const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev);
	inline double func(const double _x);
	long double Time();
	double GetSolutionX() { return solutionX; };
	double GetSolutionZ() { return solutionZ; };
};

#endif


