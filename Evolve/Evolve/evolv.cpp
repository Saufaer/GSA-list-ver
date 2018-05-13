#include "evolv.h"
void evolv::transform_P_to_D()
{
	//if (N == 1) return;
	// transformation from hypercube P to hyperinterval D
	for (int i = 0; i < N; i++)
		y[i] = y[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2;
}

void evolv::CalculateNode(double is, int n, int *u, int *v, int *l)
// вычисление вспомогательного центра u(s) и соответствующих ему v(s) и l(s)
// calculate u=u[s], v=v[s], l=l[s] by is=s
{
	int n1, i, j, k1, k2, iq;
	double iff;
	double nexp;

	iq = 1;
	for (nexp = 1, i = 0; i < n; nexp += nexp, i++);
	n1 = n - 1;
	*l = 0;
	if (is == 0)
	{
		*l = n1;
		for (i = 0; i < n; i++)
		{
			u[i] = -1;
			v[i] = -1;
		}
	}
	else if (is == (nexpExtended - extOne))
	{
		*l = n1;
		u[0] = 1;
		v[0] = 1;
		for (i = 1; i < n; i++)
		{
			u[i] = -1;
			v[i] = -1;
		}
		v[n1] = 1;
	}
	else
	{
		iff = nexpExtended;
		k1 = -1;
		for (i = 0; i < n; i++)
		{
			iff = iff / 2;
			if (is >= iff)
			{
				if ((is == iff) && (is != extOne))
				{
					*l = i;
					iq = -1;
				}
				is -= iff;
				k2 = 1;
			}
			else
			{
				k2 = -1;
				if ((is == (iff - extOne)) && (is != extNull))
				{
					*l = i;
					iq = 1;
				}
			}
			j = -k1 * k2;
			v[i] = j;
			u[i] = j;
			k1 = k2;
		}
		v[*l] = v[*l] * iq;
		v[n1] = -1 * v[n1];
	}
}

double* evolv::GetYOnX(const double _x)
{
	if (N == 1)
	{
		y[0] = _x - 0.5;
		return y;
	}

	int iu[MaxDim];
	int iv[MaxDim];
	int l;
	double d;
	int mn;
	double r;
	int iw[MaxDim];
	int it, i, j;
	double is;

	d = _x;
	r = 0.5;
	it = 0;
	mn = m * N;
	for (i = 0; i < N; i++)
	{
		iw[i] = 1;
		y[i] = 0.0;
	}
	for (j = 0; j < m; j++)
	{
		if (_x == extOne)
		{
			is = nexpExtended - extOne;
			d = extNull; // d = 0.0;
		}
		else
		{
			//Код из старой версии - уточнить работоспособность при N > 32
			d *= nexpExtended;
			is = (int)d;
			d -= is;
		}
		CalculateNode(is, N, iu, iv, &l);
		i = iu[0];
		iu[0] = iu[it];
		iu[it] = i;
		i = iv[0];
		iv[0] = iv[it];
		iv[it] = i;
		if (l == 0)
			l = it;
		else if (l == it)
			l = 0;
		r *= 0.5;
		it = l;
		for (i = 0; i < N; i++)
		{
			iu[i] *= iw[i];
			iw[i] *= -iv[i];
			y[i] += r * iu[i];
		}
	}
	return y;
}