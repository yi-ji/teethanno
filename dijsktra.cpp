#include "mex.h"
#include "math.h"
#include "matrix.h"
#include <stdlib.h>

extern const int max_len;

class Heap
{
public:
	int size;
	int* array;
	double* d;

	Heap(int capacity, double* d)
	{
		this->size = 0;
		this->array = (int*)mxMalloc(capacity*sizeof(int));
		this->d = d;
	}

	void shift_up(int pos, int node)
	{
		//mexPrintf("pos:%d node:%d\n",pos, node);
		while (pos > 0)
		{
			int parent = (pos-1)/2;
			if (this->d[this->array[parent]-1] <= this->d[node-1])
			{
				break;
			}
			this->array[pos] = this->array[parent];
			pos = parent;
		}
		this->array[pos] = node;
		//mexPrintf("pos:%d node:%d\n",pos, node);
	}

	void shift_down(int pos, int node)
	{
		int half = this->size/2;
		while (pos < half)
		{
			int child = 2*pos+1;
			int right = child + 1;
			if (right <= this->size && this->d[this->array[right]-1] < this->d[this->array[child]-1])
			{
				child = right;
			}
			if (this->d[node-1] <= this->d[this->array[child]-1])
			{
				break;
			}
			this->array[pos] = this->array[child];
			pos = child;
		}
		this->array[pos] = node;
	}

	void insert(int node)
	{
		this->shift_up(this->size,node);
		this->size += 1;
	}

	void adjust(int node)
	{
		for (int i = 0; i < this->size;i++)
		{
			if (this->array[i] == node)
			{
				this->shift_up(i,node);
				break;
			}
		}
	}

	int pop()
	{
		int top = this->array[0];
		this->size -= 1;
		this->shift_down(0,this->array[this->size]);
		return top;
	}
};

double* getE(double* E, int x, int y)
{
	return E + (x-1)*max_len*2 + y;
}

void print_heap(Heap* heap)
{
	for (int i = 0;i < heap->size;i++)
	{
		mexPrintf("%d ",heap->array[i]);
	}
	mexPrintf("end\n");
}

void init(int* done, int V_num, double* d, int* pre, int s, double* E, int* E_size, Heap* heap)
{
	for (int i = 0;i < V_num;i++)
	{
		d[i] = -1;
		pre[i] = -1;
	}
	d[s-1] = 0;

	//mexPrintf("E_size(s):%d\n",E_size[s]);

	for (int i = 0;i < E_size[s-1];i++)
	{
		int next = (int) *getE(E,s,i*2);
		//mexPrintf("next:%d\n",next);
		d[next-1] = *getE(E,s,i*2+1);
		heap->insert(next);
		*done += 1;
	}
}

void dijsktra(int* done, int V_num, double* d, int* pre, int s, double* E, int* E_size, Heap* heap)
{
	while (*done < V_num && heap->size > 0)
	{
		int node = heap->pop();
		//mexPrintf("done:%d size:%d node:%d\n", *done, heap->size, node);
		//print_heap(heap);
		*done += 1;
		for (int i = 0;i < E_size[node-1];i++)
		{
			int next = ((int) *getE(E,node,i*2));
			double dist = *getE(E,node,i*2+1);
			//mexPrintf("node:%d child:%d\n",node,next);
			if (d[next-1] < 0)
			{
				d[next-1] = d[node-1] + dist;
				pre[next-1] = node;
				heap->insert(next);
			}
			else if (d[next-1] > d[node-1] + dist)
			{
				d[next-1] = d[node-1] + dist;
				pre[next-1] = node;
				heap->adjust(next);
			}
		}
	}
}

void print(int V_num, double* d, Heap* heap)
{
	int counter = 0;
	for (int i = 0;i < V_num;i++)
	{
		if (d[i] >= 0)
		{
			mexPrintf("%f ",d[i]);
			counter += 1;
		}
	}
	mexPrintf("counter:%d\n",counter);
	mexPrintf("\n-----\n");
	mexPrintf("heap->size:%d\n",heap->size);
	print_heap(heap);
}

int trace_len(int* pre, int t)
{
	int res = 0,temp = t;
	while (temp != -1)
	{
		temp = pre[temp-1];
		res += 1;
	}
	return res;
}

void trace(int* traced_nodes, int* V, int* pre, int t)
{
	int temp = t,i = 0;
	while (temp != -1)
	{
		//mexPrintf("%d %d\n",temp-1,V[temp-1]);
		traced_nodes[i] = temp;
		//traced_nodes[i] = V[temp-1];
		i += 1;
		temp = pre[temp-1];
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
	if (nrhs == 4)
	{
		int V_num = mxGetScalar(prhs[0]);
		double* E = mxGetPr(prhs[1]);
		int* E_size = (int*) mxGetPr(prhs[2]);
		int s = mxGetScalar(prhs[3]);

		plhs[0] = mxCreateDoubleMatrix(1,V_num,mxREAL);
		double* d = mxGetPr(plhs[0]);

		plhs[1] = mxCreateNumericMatrix(1,V_num,mxINT32_CLASS,mxREAL);
		int* pre = (int*)mxGetData(plhs[1]);

		Heap heap(V_num,d);
		//mexPrintf("step 1 done\n");

		int done = 0;

		init(&done,V_num,d,pre,s,E,E_size,&heap);
		//mexPrintf("step 2 done\n");
		//print(V_num,d,&heap);
		dijsktra(&done,V_num,d,pre,s,E,E_size,&heap);
	}
	else if (nrhs == 3)
	{
		int* V = (int*) mxGetPr(prhs[0]);
		int* pre = (int*) mxGetPr(prhs[1]);
		int t = mxGetScalar(prhs[2]);

		int len = trace_len(pre,t);

		plhs[0] = mxCreateNumericMatrix(1,len,mxINT32_CLASS,mxREAL);
		int* traced_nodes = (int*)mxGetData(plhs[0]);

		trace(traced_nodes,V,pre,t);
	}
}
