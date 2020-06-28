#include <cstdio>
#include <math.h>
#include <chrono>
using namespace std;

//initialize
const int nx = 41;
const int ny = 41;
const int nt = 500;
const int nit = 50;
const int length = nx * ny;
//const int c = 1;
const int M = 1024;
const int N = (length + M - 1) / M;

__global__ void initialize_block(double *u, double *v, double *p, double *b) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= length)
		return;
	u[index] = 0;
	v[index] = 0;
	p[index] = 0;
	b[index] = 0;

}


//brackets
__global__ void build_up_b_block(double *u, double *v, double *b, double dx, double dy, double rho,double dt) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	//if (i >= length)
		//return;
	if (index < nx)
		return;
	if (index >= nx * (ny - 1))
		return;
	if (index%nx == 0)
		return;
	if (index%nx == nx - 1)
		return;
	b[index]= (rho * (1 / dt *
		((u[index + 1] - u[index - 1]) /
		(2 * dx) + (v[index + nx] - v[index - nx]) / (2 * dy)) -
		pow((u[index + 1] - u[index - 1]) / (2 * dx), 2) -
		2 * ((u[index + nx] - u[index - nx]) / (2 * dy) *
		(v[index + 1] - v[index - 1]) / (2 * dx)) -
		pow((v[index + nx] - v[index - nx]) / (2 * dy), 2)));
}


__global__ void pressure_poisson_block(double *p,double *pn, double *b, double dx, double dy, double rho) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index < nx)
		return;
	if (index >= nx * (ny - 1))
		return;
	if (index%nx == 0)
		return;
	if (index%nx == nx - 1)
		return;
	p[index] = (((pn[index + 1] + pn[index - 1]) * pow(dy, 2) +
		(pn[index + nx] + pn[index - nx]) * pow(dx, 2)) /
		(2 * (pow(dx, 2) + pow(dy, 2))) -
		pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
		b[index]);
}

//copy the array
 __global__ void copy_block(double *a, double *an) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >=length)
		return;
	an[index] = a[index];
}

 __global__ void boundary_p_block(double *p) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index >= length)
		return;
	if (index%nx == nx - 1)
		p[index] = p[index - 1];
	if (index%nx == 0)
		p[index] = p[index + 1];
	if (index < nx)
		p[index] = p[index + nx];
	if (index >= nx * (ny - 1))
		p[index] = 0;
		
}
void pressure_poisson(double *p, double dx, double dy, double *b,double *pn,double rho) {
	for (int k = 0; k < nit; k++)
	{
		copy_block<<<N,M>>>(p, pn);
		pressure_poisson_block<<<N,M>>>(p, pn, b, dx, dy, rho);
		boundary_p_block<<<N,M>>>(p);
	}
}

 __global__ void cavity_flow_block(double *u, double *v, double dt, double dx, double dy,
	double *p, double rho, double nu,double *un,double *vn) {
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index < nx)
		return;
	if (index >= nx * (ny - 1))
		return;
	if (index%nx == 0)
		return;
	if (index%nx == nx - 1)
		return;
	u[index] = (un[index] -
		un[index] * dt / dx *
		(un[index] - un[index - 1]) -
		vn[index] * dt / dy *
		(un[index] - un[index - nx]) -
		dt / (2 * rho * dx) * (p[index + 1] - p[index - 1]) +
		nu * (dt / pow(dx, 2) *
		(un[index + 1] - 2 * un[index] + un[index - 1]) +
			dt / pow(dy, 2) *
			(un[index + nx] - 2 * un[index] + un[index - nx])));

	v[index] = (vn[index] -
		un[index] * dt / dx *
		(vn[index] - vn[index - 1]) -
		vn[index] * dt / dy *
		(vn[index] - vn[index - nx]) -
		dt / (2 * rho * dy) * (p[index + nx] - p[index - nx]) +
		nu * (dt / pow(dx, 2) *
		(vn[index + 1] - 2 * vn[index] + vn[index - 1]) +
			dt / pow(dy, 2) *
			(vn[index + nx] - 2 * vn[index] + vn[index - nx])));
}

 __global__ void boundary_c1_block(double *p) {
	 int index = blockIdx.x * blockDim.x + threadIdx.x;
	 if (index < nx)
		 p[index] = 0;
 }
 __global__ void boundary_c2_block(double *p) {
	 int index = blockIdx.x * blockDim.x + threadIdx.x;
	 if (index >= length)
		 return;
	 if (index >= nx * (ny - 1))
		 p[index] = 0;
 }
 __global__ void boundary_c3_block(double *p) {
	 int index = blockIdx.x * blockDim.x + threadIdx.x;
	 if (index >= length)
		 return;
	 if (index%nx == 0)
		 p[index] = 0;
	 if (index%nx == nx - 1)
		 p[index] = 0;
 }
 __global__ void boundary_c4_block(double *p) {
	 int index = blockIdx.x * blockDim.x + threadIdx.x;
	 if (index >= length)
		 return;
	 if (index >= nx * (ny - 1))
		 p[index] = 1;
 }
void cavity_flow(double *u, double *v, double dt, double dx, double dy,
	double *p, double rho, double nu, double *b) {
	double *un = (double*)malloc(length * sizeof(double));
	double *vn = (double*)malloc(length * sizeof(double));
	double *pn = (double*)malloc(length * sizeof(double));
	cudaMallocManaged(&u, length * sizeof(double));
	cudaMallocManaged(&v, length * sizeof(double));
	cudaMallocManaged(&p, length * sizeof(double));
	cudaMallocManaged(&b, length * sizeof(double));
	cudaMallocManaged(&un, length * sizeof(double));
	cudaMallocManaged(&vn, length * sizeof(double));
	cudaMallocManaged(&pn, length * sizeof(double));
	initialize_block<<<N,M>>>(u, v, p, b);
	for (int k = 0; k < nt; k++)
	{
		copy_block <<<N, M >>> (u,un);
		copy_block <<<N, M >>> (v,vn);
		//build_up_b
		build_up_b_block<<<N,M>>>(u, v, b, dx, dy, rho, dt);
		pressure_poisson(p, dx, dy, b, pn,rho);
		cavity_flow_block<<<N,M>>>(u, v, dt, dx, dy, p, rho, nu, un, vn);
		
		
		boundary_c1_block <<<N, M >>> (u);
		boundary_c1_block <<<N, M >>> (v);
		boundary_c2_block <<<N, M >>> (v);
		boundary_c3_block <<<N, M >>> (u);
		boundary_c3_block <<<N, M >>> (v);
		boundary_c4_block <<<N, M >>> (u);
		cudaDeviceSynchronize();
	}
		for (int i = 0; i < ny; i++)
		{
			for (int j = 0; j < nx; j++)
			{
				int in = i * nx + j;
				printf("%f ", u[in]);
			}
			printf("\n");
		}
		cudaFree(u);
		cudaFree(v);
		cudaFree(p);
		cudaFree(b);
		cudaFree(un);
		cudaFree(vn);
		cudaFree(pn);
}


int main() {
	double dx = 2.0 / (nx - 1);
	double dy = 2.0 / (ny - 1);
	double rho = 1;
	double nu = 0.1;
	double dt = 0.001;
	double u[length];
	double v[length];
	double p[length];
	double b[length];
	auto tic = chrono::steady_clock::now();
	cavity_flow(u, v, dt, dx, dy, p, rho, nu, b);
	auto toc = chrono::steady_clock::now();
	double time = chrono::duration<double>(toc - tic).count();
	printf("===========\n");
	printf("above is matrix u\n");
	printf("===========\n");
	printf("runtime: %f", time);
	printf("\n===========");
}
