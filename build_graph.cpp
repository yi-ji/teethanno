#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

extern const double angle_bound;
extern const int max_len;

double* getV(double* V, int i);
int* getV(int* V, int i);
double dist(double* a, double* b);
double norm(double* a);
double dot_product(double* a, double* b);
double cal_angle(double* point, double* center, double center_dist, double* std_vec, double norm_std_vec, double angle_bound);
double fnval(const mxArray* pp, double angle);
int d2i(double a);

double* getE(double* E, int x, int y)
{
	return E + x*max_len*2 + y;
}

void build_each(int a, int b, double* center_coef, double* curv, double* graph, int* graph_size, int switcher)
{
	a -= 1;b -= 1;
	if (graph_size[a] == max_len) {return;}
	for (int i = 0;i < graph_size[a];i++)
	{
		double* p = getE(graph,a,i*2);
		if (fabs(*p-double(b+1)) < 0.001) {return;}
	}
	//mexPrintf("%d %d\n",a,b);
	double* p = getE(graph,a,graph_size[a]*2);
	*p = b+1;
	switch (switcher)
	{
		case 4:
			*(p+1) = fabs(center_coef[a]+center_coef[b]);
			break;
		case 3:
			*(p+1) = fabs((center_coef[a]+center_coef[b])/sqrt(curv[a]*curv[b]));
			break;
		default:
			*(p+1) = fabs((center_coef[a]+center_coef[b])/(curv[a]*curv[b]));
			break;
	}
	graph_size[a] += 1;
}

void build_graph(int* F, int F_num, int V_num, double* curv, double* center_dist, double* border, double* graph, int* graph_size, int switcher)
{
	//double center_basis = 10.0;
	//double* center_bias = (double*) mxMalloc(2*sizeof(double));
	//center_bias[0] = 15.2/*16*/;center_bias[1] = 23.0;
	double bias = -2.0;
	double* center_coef = (double*) mxMalloc(V_num*sizeof(double));

	switch (switcher)
	{
		case 1:
			for (int i = 0;i < V_num;i++)
			{
				if (center_dist[i] > border[i]-bias)
				{
					center_coef[i] = 10000;
				}
				else
				{
					center_coef[i] = 1;//pow((border[i] - center_dist[i]) / border[i], 4);
				}
			}
			break;
		case 2:
			for (int i = 0;i < V_num;i++)
			{
				if (center_dist[i] < border[i])
				{
					center_coef[i] = 100;
				}
				else
				{
					center_coef[i] = 1;//pow((center_dist[i] - border) / border, 0.5);
				}
			}		
			break;
		case 3:
			for (int i = 0;i < V_num;i++)
			{
				center_coef[i] = 1;
			}
			break;
		case 4:
			for (int i = 0;i < V_num;i++)
			{
				center_coef[i] = 1;
			}
			for (int i = 0;i < d2i(border[0]);i++)
			{
				center_coef[d2i(center_dist[i])-1] = 9999;
			}
			break;
		default: break;
	}

	for (int i = 0;i < F_num;i++)
	{
		int* temp = getV(F,i);
		//mexPrintf("%d\n",i);
		build_each(*temp, *(temp+1), center_coef, curv, graph, graph_size, switcher);
		build_each(*(temp+1), *temp, center_coef, curv, graph, graph_size, switcher);
		build_each(*temp, *(temp+2), center_coef, curv, graph, graph_size, switcher);
		build_each(*(temp+2), *temp, center_coef, curv, graph, graph_size, switcher);
		build_each(*(temp+1), *(temp+2), center_coef, curv, graph, graph_size, switcher);
		build_each(*(temp+2), *(temp+1), center_coef, curv, graph, graph_size, switcher);
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	int* F = (int*) mxGetPr(prhs[0]);
	int F_num = mxGetN(prhs[0]);
	double* curv = mxGetPr(prhs[1]);
	double* center_dist = mxGetPr(prhs[2]);
	int V_num = mxGetM(prhs[1]);
	double* border = mxGetPr(prhs[3]);
	int switcher = mxGetScalar(prhs[4]);

	plhs[0] = mxCreateDoubleMatrix(max_len*2,V_num,mxREAL);
	double* graph = mxGetPr(plhs[0]);
	plhs[1] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
	int* graph_size = (int*)mxGetData(plhs[1]);

	build_graph(F,F_num,V_num,curv,center_dist,border,graph,graph_size,switcher);
}








