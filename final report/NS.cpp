#include <iostream>
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

void initialize(double *u, double *v, double *p, double *b) {
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			u[i*nx + j] = 0;
			v[i*nx + j] = 0;
			p[i*nx + j] = 0;
			b[i*nx + j] = 0;
		}
	}
}

void build_up_b(double *b, double rho, double dt, double *u, double *v, double dx, double dy) {
	for (int i = 1; i < ny-1; i++)
	{
		for (int j = 1; j < nx-1; j++)
		{
			int index = nx * i + j;
			b[index] = (rho * (1 / dt *
				((u[index + 1] - u[index - 1]) /
				(2 * dx) + (v[index + nx] - v[index - nx]) / (2 * dy)) -
				pow((u[index + 1] - u[index - 1]) / (2 * dx), 2) -
				2 * ((u[index + nx] - u[index - nx]) / (2 * dy) *
				(v[index + 1] - v[index - 1]) / (2 * dx)) -
				pow((v[index + nx] - v[index - nx]) / (2 * dy), 2)));
		}
	}
}

void pressure_poisson(double *p, double dx, double dy, double *b) {
	double *pn = (double*)malloc(length * sizeof(double));
	for (int k = 0; k < nit; k++)
	{
		for (int i = 0; i < length; i++)
			pn[i] = p[i];
		for (int i = 1; i < ny-1; i++)
		{
			for (int j = 1; j < nx-1; j++)
			{
				int index = nx * i + j;
				p[index] = (((pn[index + 1] + pn[index - 1]) * pow(dy, 2) +
					(pn[index + nx] + pn[index - nx]) * pow(dx, 2)) /
					(2 * (pow(dx, 2) + pow(dy, 2))) -
					pow(dx, 2) * pow(dy, 2) / (2 * (pow(dx, 2) + pow(dy, 2))) *
					b[index]);
			}
		}
		// dp/dx = 0 at x = 2 and x = 0
		for (int i = 0; i < ny; i++)
			p[(i + 1)*nx - 1] = p[(i + 1)*nx - 2];
		for (int i = 0; i < ny; i++)
			p[i*nx] = p[i*nx + 1];
		// dp/dy = 0 at y = 0
		for (int i = 0; i < nx; i++)
			p[i] = p[i + nx];
		for (int i = nx*(ny-1); i < nx*ny; i++)
			p[i] = 0;
		

	}
}
void cavity_flow( double *u, double *v, double dt, double dx, double dy,
	double *p, double rho, double nu,double *b) {
	int length = nx * ny;
	double *un = (double*)malloc(length * sizeof(double));
	double *vn = (double*)malloc(length * sizeof(double));
	for (int k = 0; k < nt; k++)
	{
		for (int i = 0; i < length; i++)
		{
			un[i] = u[i];
			vn[i] = v[i];
		}
		build_up_b(b, rho, dt, u, v, dx, dy);
		pressure_poisson(p, dx, dy, b);
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				int index = i * nx + j;
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
		}
		for (int i = 0; i < nx; i++)
		{
			u[i] = 0;
			v[i] = 0;
			v[nx*(ny - 1) + i] = 0;
		}

		for (int i = 0; i < ny; i++)
		{
			u[i*nx] = 0;
			v[i*nx] = 0;
			u[(i + 1)*nx - 1] = 0;
			v[(i + 1)*nx - 1] = 0;
		}

		for (int i = 0; i < nx; i++)
			u[nx*(ny - 1) + i] = 1;
	}
}


int main() {
	double dx = 2.0/(nx-1);
	double dy = 2.0/(ny-1);
	double rho = 1;
	double nu = 0.1;
	double dt = 0.001;
	double u[length];
	double v[length];
	double p[length];
	double b[length];
	initialize(u, v, p, b);
	auto tic = chrono::steady_clock::now();
	cavity_flow(u, v, dt, dx, dy, p, rho, nu, b); 
	auto toc = chrono::steady_clock::now();
	double time = chrono::duration<double>(toc - tic).count();
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			int in = i * nx + j;
			printf("%f ", u[in]);
		}
		printf("\n");
	}
	printf("===========\n");
	printf("above is matrix u\n");
	printf("===========\n");
	printf("runtime: %f", time);
	printf("\n===========");
}