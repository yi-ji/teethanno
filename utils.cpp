#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

extern const double angle_bound = 1.85;
extern const int max_len = 20;

double* getV(double* V, int i)
{
	return V + i*3;
}

int* getV(int* V, int i)
{
	return V + i*3;
}

double* a_to_b(double* a, double* b)
{
	double* vec = (double*) mxMalloc(3*sizeof(double));
	vec[0] = b[0] - a[0];
	vec[1] = b[1] - a[1];
	vec[2] = b[2] - a[2];
	return vec;
}

double dist(double* a, double* b)
{
	double dx = a[0] - b[0];
	double dy = a[1] - b[1];
	double dz = a[2] - b[2];
	return sqrt(dx*dx + dy*dy + dz*dz);
}

double norm(double* a)
{
	if (a == NULL) {return 0;}
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

double dot_product(double* a, double* b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double cal_cos(double* a, double* b)
{
	return dot_product(a,b)/(norm(a)*norm(b));
}

double cal_cos_xy(double* a, double* b)
{
	double* c = (double*) mxMalloc(3*sizeof(double));
	double* d = (double*) mxMalloc(3*sizeof(double));
	c[0] = a[0];c[1] = 0.0;c[2] = a[2];
	d[0] = b[0];d[1] = 0.0;d[2] = b[2];
	return cal_cos(c,d);
}

double cal_angle(double* point, double* center, double center_dist, double* std_vec, double norm_std_vec, double angle_bound)
{
	double* vec = a_to_b(center,point);
	double angle = acos(dot_product(vec,std_vec)/(center_dist*norm_std_vec));
	angle = fmax(fmin(angle,angle_bound),-angle_bound);
	return angle;
}

double fnval(const mxArray* pp, double angle)
{
	mxArray* ans[1];
	mxArray* input[2];
	input[0] = (mxArray*) pp;
	input[1] = mxCreateDoubleScalar(angle);
	mexCallMATLAB(1,ans,2,input,"fnval");
	return mxGetScalar(ans[0]);
}

int* double_capacity(int* a, int len)
{
	int* b = (int*) mxMalloc(len*2*sizeof(int));
	for (int i = 0;i < len;i++)
	{
		b[i] = a[i];
	}
	return b;
}

int d2i(double a)
{
	int b = int(a);
	if (fabs(double(b)-a) < 0.01)
	{
		return b;
	}
	else
	{
		return b+1;
	}
}