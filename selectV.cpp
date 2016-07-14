#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

const int max_len = 100;

int* double_capacity(int* a, int len);

int d2i(double a);

void mycopy(bool* a, bool* b, int len)
{
	for (int i = 0;i < len;i++)
	{
		b[i] = a[i];
	}
}

int visit(int node, double* E, int* E_size, bool* trimed, bool* visited)
{
	int queue_size = 100;
	int* q = (int*) mxMalloc(queue_size*sizeof(int));
	int head = 0, tail = 0;
	q[head] = node;
	head += 1;
	visited[node] = true;
	trimed[node] = true;

	while (tail < head)
	{
		int temp = q[tail];
		for (int i = 0;i < E_size[temp];i++)
		{
			int child = d2i(E[temp*max_len+i])-1;
			//mexPrintf("temp:%d child:%d\n",temp,child);
			if (visited[child]) {continue;}
			visited[child] = true;
			trimed[child] = true;
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

	return tail+1;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{

	double* E = mxGetPr(prhs[0]);
	int* E_size = (int*) mxGetPr(prhs[1]);
	double* nodes = mxGetPr(prhs[2]);
	int nodes_num = mxGetScalar(prhs[3]);
	int step_threshold = mxGetScalar(prhs[4]);

	bool* trimed = (bool*) mxMalloc(nodes_num*sizeof(bool));
	bool* trimed_bak = (bool*) mxMalloc(nodes_num*sizeof(bool));
	bool* visited = (bool*) mxMalloc(nodes_num*sizeof(bool));

	for (int i = 0;i < nodes_num;i++)
	{
		trimed[i] = false;
		trimed_bak[i] = false;
		visited[i] = false;
	}

	for (int i = 0;i < nodes_num;i++)
	{
		if (visited[i]) {continue;}
		//mexPrintf("%d begin\n",i);
		int size = visit(i,E,E_size,trimed,visited);
		//if (size > 2) {mexPrintf("size:%d\n",size);}
		//mexPrintf("%d done\n",i);
		if (size >= step_threshold)
		{
			mycopy(trimed,trimed_bak,nodes_num);
		}
		else
		{
			mycopy(trimed_bak,trimed,nodes_num);
		}
	}

	int V_num = 0;

	for (int i = 0;i < nodes_num;i++)
	{
		if (trimed[i]) {V_num += 1;}
	}

	plhs[0] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
	int* V = (int*)mxGetData(plhs[0]);

	int j = 0;
	for (int i = 0;i < nodes_num;i++)
	{
		if (trimed[i])
		{
			V[j] = d2i(nodes[i]);
			j += 1;
		}
	}

	plhs[1] = mxCreateDoubleScalar(double(V_num));

}