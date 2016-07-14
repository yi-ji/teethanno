#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

extern const int max_len;

int* double_capacity(int* a, int len);

void dye(int* graph, int* graph_size, bool* flag, int* result, int center, int color)
{
	int queue_size = 100;
	int* q = (int*) mxMalloc(queue_size*sizeof(int));
	int head = 0, tail = 0;
	q[head] = center - 1;
	head += 1;
	result[center-1] = color;

	while (tail < head)
	{
		int temp = q[tail];
		for (int i = 0;i < graph_size[temp];i++)
		{
			int child = graph[(temp*max_len+i)*2] - 1;
			//mexPrintf("temp:%d child:%d tail:%d head:%d size:%d\n",temp,child,tail,head,queue_size);
			result[child] = color;
			if (flag[child]) {continue;}
			flag[child] = true;
			if (head >= queue_size)
			{
				q = double_capacity(q,queue_size);
				queue_size *= 2;
			}
			q[head] = child;
			head += 1;
		}
		tail += 1;
	}

	mxFree(q);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	int* graph = (int*) mxGetPr(prhs[0]);
	int* graph_size = (int*) mxGetPr(prhs[1]);
	int V_num = mxGetScalar(prhs[2]);
	int* boundary = (int*) mxGetPr(prhs[3]);
	int boundary_len = mxGetN(prhs[3]);
	int center = mxGetScalar(prhs[4]);
	int color = mxGetScalar(prhs[5]);

	bool* flag = (bool*) mxMalloc(V_num*sizeof(bool));

	plhs[0] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
	int* result = (int*)mxGetData(plhs[0]);

	for (int i = 0;i < V_num;i++)
	{
		flag[i] = false;
		result[i] = 0;
	}

	for (int i = 0;i < boundary_len;i++)
	{
		flag[boundary[i]-1] = true;
	}

	dye(graph,graph_size,flag,result,center,color);
}