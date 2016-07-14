#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

extern const double angle_bound;
const int max_len = 100;

double* getV(double* V, int i);
double dist(double* a, double* b);
double norm(double* a);
double dot_product(double* a, double* b);
double cal_angle(double* point, double* center, double center_dist, double* std_vec, double norm_std_vec, double angle_bound);
double fnval(const mxArray* pp, double angle);

double* getE(double* E, int x, int y)
{
	return E + x*max_len*2 + y;
}

void buildE(double* E, int* E_size, int V_num, double* V, double* curv, const mxArray* inner_pp, const mxArray* outer_pp, double dist_threshold, double km_exp, double dist_exp, double weight_threshold, int switcher, double* center, double* std_vec, int s, int t)
{
	double center_basis = 10.0, bad_case_range = 7.0;
	double* center_bias = (double*) mxMalloc(2*sizeof(double));
	center_bias[0] = 16.0;center_bias[1] = 23.0; // need to discard
	double norm_std_vec = norm(std_vec);

	//mexPrintf("step 3 done\n");

	for (int i = 0;i < V_num;i++)
	{
		double center_coef = 1.0,center_dist,angle,dist1,dist2;
		switch (switcher)
		{
			case 0: break;
			case 1:
				center_coef = (dist(center,getV(V,i))-center_bias[switcher-1]) / center_basis;
				center_coef = pow(center_coef,2);
				//if (dist(getV(V,t),getV(V,i)) < bad_case_range || dist(getV(V,i),getV(V,s)) < bad_case_range) {center_coef = pow(center_coef,7);}
				break;
			case 2:
				center_coef = center_basis / (dist(center,getV(V,i))-center_bias[switcher-1]);
				center_coef = pow(center_coef,2);
				//if (dist(getV(V,t),getV(V,i)) < bad_case_range || dist(getV(V,i),getV(V,s)) < bad_case_range) {center_coef = pow(center_coef,7);}
				break;
			case 3:	
				center_dist = dist(center,getV(V,i));
				angle = cal_angle(getV(V,i),center,center_dist,std_vec,norm_std_vec,angle_bound);
				dist1 = fnval(inner_pp,angle);
				dist2 = fnval(outer_pp,angle);
				if (center_dist > (dist1*2.0+dist2)/3.0)
				{
					center_coef = 999.0;
				}
				else
				{
					center_coef = pow(((dist1+dist2)/2.0 - center_dist) / center_dist, 4);
				}
				break;
			case 4:
				center_dist = dist(center,getV(V,i));
				angle = cal_angle(getV(V,i),center,center_dist,std_vec,norm_std_vec,angle_bound);
				dist1 = fnval(inner_pp,angle);
				dist2 = fnval(outer_pp,angle);
				if (center_dist < (dist1+dist2*2.0)/3.0)
				{
					center_coef = 999.0;
				}
				else
				{
					center_coef = pow((center_dist - (dist1+dist2)/2.0) / center_dist, 4);
					//center_coef = line_dist < 1.5 ? center_coef*line_dist*0.02 : center_coef;
					if (dist(getV(V,i),getV(V,t)) < bad_case_range || dist(getV(V,i),getV(V,s)) < bad_case_range) {center_coef = pow(1.0/center_coef,2);}
				}
				break;
			case 5:
				break;
			default: break;
		}
		for (int j = 0;j < V_num;j++)
		{
			if (i == j || E_size[i] >= max_len) {continue;}
			double temp = dist(getV(V,i),getV(V,j));
			if (temp >= dist_threshold) {continue;}
			//mexPrintf("step 4 done\n");

			//if (switcher == 3 || switcher == 4) {center_coef *= temp > 0.5 ? (temp-0.5)*100 : 1.0;}

			if (switcher == 0)
			{
				//mexPrintf("%d %d\n",i,j);
				double* p = E + max_len*i + E_size[i];
				*p = j+1;
				E_size[i] += 1;
				continue;
			}

			temp = pow(temp,dist_exp) / pow(fabs(curv[i]),km_exp);
			//mexPrintf("step 4.5 done\n");
			if (temp >= weight_threshold) {continue;}
			
			//mexPrintf("step 5 done\n");
			
			//mexPrintf("step 6 done\n");
			//mexPrintf("%d %d %f\n",i,E_size[i],temp*center_coef);
			double* p = getE(E,i,E_size[i]*2);
			*p = j+1; *(p+1) = fabs(temp*center_coef);
			E_size[i] += 1;
			//mexPrintf("step 7 done\n");
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if (nrhs > 10)
	{
		int V_num = mxGetScalar(prhs[0]);
		double* V = mxGetPr(prhs[1]);
		double* curv = mxGetPr(prhs[2]);
		const mxArray* inner_pp = prhs[3];
		const mxArray* outer_pp = prhs[4];
		double dist_threshold = mxGetScalar(prhs[5]);
		double km_exp = mxGetScalar(prhs[6]);
		double dist_exp = mxGetScalar(prhs[7]);
		double weight_threshold = mxGetScalar(prhs[8]);
		int switcher = mxGetScalar(prhs[9]);
		double* center = mxGetPr(prhs[10]);
		double* std_vec = mxGetPr(prhs[11]);
		int s = mxGetScalar(prhs[12]) - 1;
		int t = mxGetScalar(prhs[13]) - 1;

	//mexPrintf("step 1 done\n");

		plhs[0] = mxCreateDoubleMatrix(max_len*2,V_num,mxREAL);
		double* E = mxGetPr(plhs[0]);

		plhs[1] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
		int* E_size = (int*)mxGetData(plhs[1]);

		//mexPrintf("step 2 done\n");
		//mexPrintf("s:%f %f %f\n",getV(V,s)[0],getV(V,s)[1],getV(V,s)[2]);
		//mexPrintf("t:%f %f %f\n",getV(V,t)[0],getV(V,t)[1],getV(V,t)[2]);
		buildE(E, E_size, V_num, V, curv, inner_pp, outer_pp, dist_threshold, km_exp, dist_exp, weight_threshold, switcher, center, std_vec, s, t);
	}
	else
	{
		int V_num = mxGetScalar(prhs[0]);
		double* V = mxGetPr(prhs[1]);
		double dist_threshold = mxGetScalar(prhs[2]);

		plhs[0] = mxCreateDoubleMatrix(max_len,V_num,mxREAL);
		double* E = mxGetPr(plhs[0]);

		plhs[1] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
		int* E_size = (int*)mxGetData(plhs[1]);


		buildE(E, E_size, V_num, V, NULL, NULL, NULL, dist_threshold, 0, 0, 0, 0, NULL, NULL, 0, 0);
	}
}








