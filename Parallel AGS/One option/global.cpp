#include "global.h"

global::global()
{
	left = 0;
	right = 0;
	r = 0;
	time = 0.L;
}

Pointer global::solve(const double _left, const double _right, const double _r, const int _N_max, const double _Eps)
{
#pragma omp single

	Pointer p;


	double solutionX = 0;
	double solutionZ = 0;

	left = _left;
	double curr_left = _left;
	right = _right;
	double curr_right = _right;
	r = _r;
	
	Values_Z.clear();
	Arguments_X.clear();

	double z_begin = func(left); 
	double z_end = func(right); /
								 
	double temp = 0.; 
	int k = 1; 
	int new_size;
	double M_big, m_small = 1.; /
	int R_max_index = 0;

	
	double new_point;
						
						
	std::list<double> z;
	std::list<double> x;
	std::list<double> M_vector;
	std::list<double> R_vector;

								
	std::list<double>::iterator place_x, place_z, place_M, place_R;
	std::list<double>::iterator iter_x, iter_z, iter_R, iter_M;
	////
	left = curr_left;
	right = curr_right;
	new_point = temp = (curr_left + curr_right) / 2;

	x.push_back(left);
	x.push_back(new_point);
	x.push_back(right);

	iter_x = x.begin();
	for (int i = 0; i < 3; i++)
	{
		z.push_back(func(*iter_x++));
	}

	iter_z = z.begin();
	iter_x = x.begin();


	block val_z; 
	block arg_x; 

	{

		//поиск M
		for (int i = 1; i < 3; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;

			M_vector.push_back(M(val_z.right, val_z.left, arg_x.right, arg_x.left));
		}
		double max = *M_vector.begin();
		if (max < *M_vector.begin() + 1) max = *M_vector.begin() + 1;
		M_big = max;
		if (M_big > 0.) m_small = r * M_big;
		else if (M_big == 0.) m_small = 1.;

		R_max_index = 0;
		max = 0.;

		iter_z = z.begin();
		iter_x = x.begin();

		
		for (int i = 1; i < 3; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;

			temp = R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left);
			R_vector.push_back(temp);

			if (max < temp)
			{
				max = temp;
				R_max_index = i;
			}
		}

		iter_z = z.begin();
		iter_x = x.begin();

		for (int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}

		solutionZ = *iter_z;
		solutionX = *iter_x;
		curr_left = *--iter_x;
		curr_right = *++iter_x;

		iter_z = z.begin();
		iter_x = x.begin();
	}


	while (k < _N_max && (fabs(curr_left - curr_right) > _Eps))
	{
		new_size = 2 + k;

		place_z = z.begin();
		place_x = x.begin();
		for (int i = 0; i < R_max_index; i++)
		{
			place_z++;
			place_x++;
		}
		val_z.right = *place_z--;
		val_z.left = *place_z++;

		left = curr_left;
		right = curr_right;
		new_point = 0.5*(right + left) - (val_z.right - val_z.left) / (2 * m_small);

		z.insert(place_z, func(new_point));
		x.insert(place_x, new_point);

		iter_z = z.begin();
		iter_x = x.begin();
		place_M = M_vector.begin();
		place_R = R_vector.begin();
		for (int i = 0; i < R_max_index - 1; i++)
		{
			iter_z++;
			iter_x++;
			if (i == R_max_index - 1) break;
			place_M++;
			place_R++;
		}

	
		for (int i = 0; i < 2; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;
			if (i == 0)
			{
				M_vector.insert(place_M, M(val_z.right, val_z.left, arg_x.right, arg_x.left));
			}
			else
			{
				*place_M = M(val_z.right, val_z.left, arg_x.right, arg_x.left);
			}
		}

		iter_M = M_vector.begin();
		double max = *iter_M++;
		temp = *iter_M;
		for (int i = 0; i < new_size - 1; i++)
		{
			if (max < temp) max = temp;
		}
		M_big = max;
		if (M_big > 0.) m_small = r * M_big;
		else if (M_big == 0.) m_small = 1.;
		R_max_index = 1;
		max = 0.;

		--iter_z; --iter_z;
		--iter_x; --iter_x;

		for (int i = 0; i < 2; i++)
		{
			val_z.left = *iter_z++;
			val_z.right = *iter_z;
			arg_x.left = *iter_x++;
			arg_x.right = *iter_x;
			if (i == 0)	R_vector.insert(place_R, R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left));
			else
				*place_R = R(m_small, arg_x.right, arg_x.left, val_z.right, val_z.left);
		}

		iter_R = R_vector.begin();
		max = *iter_R;
		for (int i = 0; i < new_size - 1; i++)
		{
			temp = *iter_R++;
			if (max < temp)
			{
				max = temp;
				R_max_index = i + 1;
			}
		}

		iter_x = x.begin();
		iter_z = z.begin();
		for (int i = 0; i < R_max_index; i++)
		{
			iter_z++;
			iter_x++;
		}

		solutionZ = *iter_z;
		solutionX = *iter_x;

		p.z = solutionZ;
		p.x = solutionX;

		curr_left = *--iter_x;
		curr_right = *++iter_x;
		++k;
	}

	--new_size;
	p.steps = new_size;
	iter_x = x.begin();
	iter_z = z.begin();
	for (register int i = 0; i < new_size; i++)
	{
		Arguments_X.push_back(*iter_x++);
		Values_Z.push_back(*iter_z++);
	}
	return p;
}

inline double global::R(const double _m_small, const double _x_curr, const double _x_prev, const double _z_curr, const double _z_prev)
{
	return _m_small * (_x_curr - _x_prev) + pow(_z_curr - _z_prev, 2) / (_m_small*(_x_curr - _x_prev)) - 2 * (_z_curr + _z_prev);
}

inline double global::M(const double _z_curr, const double _z_prev, const double _x_curr, const double _x_prev)
{
	return fabs((_z_curr - _z_prev) / (_x_curr - _x_prev));
}

inline double global::func(const double _x)
{
	return sin(_x) + sin((10. * _x) / 3.);
}

double global::GetTime()
{
	return time;
}
Pointer global::ParallelSolve(int threads, const double _left, const double _right, const double _r, const int _N_max, const double _Eps)
{

	omp_set_num_threads(threads);

	int proc_max = omp_get_max_threads();
	Pointer *Z = new Pointer[proc_max];

	double start, end;
	time = 0.L;

	start = omp_get_wtime();

#pragma omp parallel for
	for (int i = 0; i < proc_max; i++)
	{

		Z[i] = solve(_left + (_right - _left)*i / proc_max, _left + (_right - _left)*(i + 1) / proc_max, _r, _N_max, _Eps);

	}

	Pointer res;
	res.z = Z[proc_max - 1].z;
	res.x = Z[proc_max - 1].x;
	end = omp_get_wtime();
	time = end - start;
	for (int i = 0; i < proc_max; i++)
	{
		if (Z[i].z < res.z)
		{
			res.z = Z[i].z;
			res.x = Z[i].x;

		}
		steps += Z[i].steps;
	}
	return res;

}