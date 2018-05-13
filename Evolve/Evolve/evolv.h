#include <stdlib.h>
#include <stdio.h>
#include <string>
const int MaxDim = 50;
class evolv
{
public:


int      m = 10;             // accuracy of decomposition of hypercube
int      N = 1;             // dimension
int mDimension = 1;
double   A[MaxDim];     // left and
double   B[MaxDim];     // right bounds of search area

double*  y = new double[N];             // y point from hypercube [-1/2, 1/2]^N


const double extNull = 0; // = 0.0; //Extended(0.0);
const double extOne = 1;  // = 1.0; //Extended(1.0);
const double extHalf = 0.5; // = 0.5; //Extended(0.5);
double nexpExtended;


void GetBounds(double *lower, double *upper)
{
	for (int i = 0; i < mDimension; i++)
	{
		lower[i] = -2.2;
		upper[i] = 1.8;
	}
}
void transform_P_to_D();

void CalculateNode(double is, int n, int *u, int *v, int *l);

double* GetYOnX(const double _x);

void GetImage( double x, double *_y, int EvolventNum)
{
	GetBounds(A,B);
	int i;
	for (nexpExtended = extOne, i = 0; i < N; nexpExtended += nexpExtended, i++)
		;
	// в одиночной развертке EvolventNum не используется
	//   введен, чтобы работал полиморфизм в множественных развертках


	//if ((x.toDouble() < 0) || (x.toDouble() > 1))
	//{
	//	throw EXCEPTION("x is out of range");
	//}
	// x ---> y
	GetYOnX(x); // it saves return value to y, so no need to call operator= again

	transform_P_to_D();

	memcpy(_y, y, N * sizeof(double));
}

};