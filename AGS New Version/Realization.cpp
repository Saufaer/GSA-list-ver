#include "Solve.h"
Pointer global::Simple(/*std::vector<ChX> &RX,*/ double left, double right, double * coeff, int Nmax, double E)
{
	std::vector<ChX> RX;
	//вставка первых двух RX
	ChX RXpoint;
	RXpoint.x = left;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	
	
	
	RXpoint.x = right;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	Pointer point;
	point.steps = 0;
	point.x = left;
	point.z = func(left, coeff);

	double maxR;
	double currentE = right - left;

	double minz;
	double minx;
	//дл€ концов отрезков
	if (func(right, coeff) < func(left, coeff))
	{
		point.x = right;
		point.z = func(right, coeff);
		minz = point.z;
		minx = point.x;
	}
	else
	{
		point.x = left;
		point.z = func(left, coeff);
		minz = point.z;
		minx = point.x;
	}

	RX[0].R = RSimple(RX, 1);

	double newX = 0;
	while ((point.steps < Nmax) && (currentE > E))
	{
		if (minz > func(newX, coeff))
		{
			minz = func(newX, coeff);
			minx = newX;
		}
		maxR = RX[0].R;
		newX = XSimple(RX, coeff, 1);
		
		//ќпределение интервала с максимальным R
		for (int i = 1; i < RX.size(); i++)
		{
			if (RX[i - 1].R > maxR)
			{
				maxR = RX[i - 1].R;
				newX = XSimple(RX, coeff, i);
			}
			//посчЄт дл€ услови€ выхода
			if (RX[i].x - RX[i - 1].x < currentE)
			{
				currentE = RX[i].x - RX[i - 1].x;
			}
		}

		//¬ставка RX
		
		RXpoint.x = newX;
		RXpoint.R = 0;
		RX.push_back(RXpoint);

		std::sort(RX.begin(), RX.end(), obj);

		
		//новые R
		for (int i = 1; i < RX.size(); i++)
		{
			RX[i - 1].R = RSimple(RX, i);
		}
		point.steps++;

	}
	point.x = minx;
	point.z = minz;
	return point;
}

Pointer global::Piavsky(/*std::vector<ChX> &RX,*/ double left, double right, double * coeff, double r, int Kmax, double E)
{
	std::vector<ChX> RX;
	//вставка первых двух RX
	ChX RXpoint;
	RXpoint.x = left;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	RXpoint.x = right;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	Pointer point;
	point.steps = 0;
	point.x = left;
	point.z = func(left, coeff);

	double maxR;
	double currentE = right - left;
	
	//дл€ концов отрезков
	if (func(right, coeff) < func(left, coeff))
	{
		point.x = right;
		point.z = func(right, coeff);	
	}
	else
	{
		point.x = left;
		point.z = func(left, coeff);
	}

	double m = CalculateM(r, RX, coeff);
	RX[0].R = RPiavsky(RX, m, 1);
	double newX = 0;
	while ((point.steps < Kmax) && (currentE > E))
	{
		maxR = RX[0].R;
		newX = X(RX, m, coeff, 1);

		//ќпределение интервала с максимальным R
		for (int i = 1; i < RX.size(); i++)
		{
			if (RX[i - 1].R > maxR)
			{
				maxR = RX[i - 1].R;
				newX = X(RX, m, coeff, i);
			}
			if (RX[i].x - RX[i - 1].x < currentE)
			{
				//посчЄт дл€ услови€ выхода
				currentE = RX[i].x - RX[i - 1].x;
			}
		}

		//¬ставка RX

		RXpoint.x = newX;
		RXpoint.R = 0;
		RX.push_back(RXpoint);
		std::sort(RX.begin(), RX.end(), obj);


		//новые R
		m = CalculateM(r, RX, coeff);
		for (int i = 1; i < RX.size(); i++)
		{
			
			RX[i - 1].R = RPiavsky(RX, m, i);
		}
		point.steps++;
		 
	}
	point.x = newX;
	point.z = func(newX, coeff);
	return point;


}

Pointer global::Strongin( double left, double right, double * coeff, double r, int Kmax, double E)
{
	
	std::vector<ChX> RX;
	//вставка первых двух RX
	ChX RXpoint;
	RXpoint.x = left;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	RXpoint.x = right;
	RXpoint.R = 0;
	RX.push_back(RXpoint);
	Pointer point;
	point.steps = 0;
	point.x = left;
	point.z = func(left, coeff);

	double maxR;
	double currentE = right - left;
	//дл€ концов отрезков
	if (func(right, coeff) < func(left, coeff))
	{
		point.x = right;
		point.z = func(right, coeff);
	}
	else
	{
		point.x = left;
		point.z = func(left, coeff);
	}

	double m = CalculateM(r, RX, coeff);
	RX[0].R = RStrongin(RX, m, 1);
	double newX = 0;
	while ((point.steps < Kmax) && (currentE > E))
	{
		maxR = RX[0].R;
		newX = X(RX, m, coeff, 1);
		
		//ќпределение интервала с максимальным R
		for (int i = 1; i < RX.size(); i++)
		{
			if (RX[i - 1].R > maxR)
			{
				maxR = RX[i - 1].R;
				newX = X(RX, m, coeff, i);
			}
			if (RX[i].x - RX[i - 1].x < currentE)
			{
				currentE = RX[i].x - RX[i - 1].x;
			}
		}

		//¬ставка RX
		
		RXpoint.x = newX;
		RXpoint.R = 0;
		RX.push_back(RXpoint);
		std::sort(RX.begin(), RX.end(), obj);

		//новые R
		m = CalculateM(r, RX, coeff);

		for (int i = 1; i < RX.size(); i++)
		{
			RX[i - 1].R = RStrongin(RX, m, i);
			
		}
		point.X.push_back(newX);
		point.steps++;

	}
	point.x = newX;
	point.z = func(newX, coeff);
	
	return point;


}


Pointer global::ParStrongin(int threads, double left, double right, double * coeff, double r, int Kmax, double E)
{

	
	
	omp_set_num_threads(threads);

	int proc_max = threads;
	Pointer *Z = new Pointer[proc_max];



#pragma omp parallel for
	for (int i = 0; i < proc_max; i++)
	{

		Z[i] = Strongin( left + (right - left) * i / proc_max, left + (right - left)*(i + 1) / proc_max, coeff, r, Kmax, E);
	}

	Pointer res;
	res.z = Z[0].z;
	res.x = Z[0].x;
	res.steps = 0;

	for (int i = 0; i < proc_max; i++)
	{
		if (Z[i].z < res.z)
		{
			res.z = Z[i].z;
			res.x = Z[i].x;

		}
		res.steps += Z[i].steps;
		for (int j = 0; j < Z[i].X.size(); j++)
		{
			res.X.push_back(Z[i].X[j]);
		}
	}

	
	return res;
}