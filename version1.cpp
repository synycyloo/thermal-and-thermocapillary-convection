#include <iostream>
#include <vector>

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <mkl.h>
#include <mkl_pardiso.h>
#include <fstream>
#include <sstream>
#include <cmath>
//#include <numbers>

/*
using::std::numbers::pi;
using::std::sin;
using::std::cos;
using::std::exp;
*/

/* PARDISO prototype. */
#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#define PARDISO pardiso_
#endif

double left_x = 0, left_y = 0; // левая граница
double right_x, right_y; // правая граница
//double tau = 0.005; // для re =1000
double tau = 0.01; // шаг по времени
double nu; // коэффициент вязкости 

// без учета du/dy = Ma dT/dx:
//double pr = 1.0; // число Прандтля 
//double gr = 1000; // число Грасгофа

double Le = 1.0; // число Льюиса
double R = 1.0;

double pr = 5.0;
double gr = 0.0;
//double gr = 5.0;
double Ma = 1000.0; // число Марангони 
double Ma_C = Ma * Le / R; // число Марангони для C
//double Ma_C = 0.0; // число Марангони для C
double Sc = pr* Le; // число Шмидта
double step_x, step_y; // шаг по x и y
double u_top = 0.0;
size_t nx = 100, ny = 50; // количество ячеек сетки 
std::vector<double> nodes_x, nodes_y; // координаты узлов сетки
std::vector<double> h_x, h_y; // шаги сетки по осям 
std::vector<double> bar_x, bar_y; // шаги сетки в точках i+1/2, j+1/2
std::vector<double> zeros(std::max(nx, ny) + 1, 0);

std::vector<double> u, v, p;

std::vector<double> T, C;

// общее число неизвестных
size_t size_m = 2 + 8 * (nx - 1) + (2 + 11 * (nx - 1)) * (ny - 2) + 2 + 8 * (nx - 1) + // u_block
nx + (10 + 10 + 11 * (nx - 2)) * (ny - 1) + nx + // v_block
1 + 4 * (nx - 1) + 4 * nx * (ny - 1);

size_t size_c = (nx + 1) * ny + nx * (ny + 1) + nx * ny;

double* c_m = new double[size_m];					// matrix
int* c_r = new int[size_c + 1];						// row
int* c_c = new int[size_m];							// col
double* c_b = new double[size_c];					// right part
double* c_s = new double[size_c];					// solution

// -----для T-----
// общее число неизвестных
size_t size_m_T = 2 + 4 * (nx - 1) + (ny - 1) * (2 + 5 * (nx - 1)) + 2 + 4 * (nx - 1);


size_t size_c_T = (nx + 1) * (ny + 1);

double* c_m_T = new double[size_m_T];					// matrix
int* c_r_T = new int[size_c_T + 1];						// row
int* c_c_T = new int[size_m_T];							// col
double* c_b_T = new double[size_c_T];					// right part
double* c_s_T = new double[size_c_T];					// solution

// -----для C-----
// общее число неизвестных
size_t size_m_C = 2 + 4 * (nx - 1) + (ny - 1) * (2 + 5 * (nx - 1)) + 2 + 4 * (nx - 1);


size_t size_c_C = (nx + 1) * (ny + 1);

double* c_m_C = new double[size_m_C];					// matrix
int* c_r_C = new int[size_c_C + 1];						// row
int* c_c_C = new int[size_m_C];							// col
double* c_b_C = new double[size_c_C];					// right part
double* c_s_C = new double[size_c_C];					// solution

std::string PATH = ".\\";

/*
void read_data() {
	std::ifstream input("input_data.txt");
	input >> right_x >> right_y >> nu;

	std::cout << right_x << "\t" << right_y << "\t" << nu << "\n";
}
*/

void read_data() {
	right_x = 1.0;
	right_y = 0.5;
	nu = 1.0;
}

// initialization

void init_u() {
	u.resize((nx + 1) * ny);
	// fill u
	for (size_t j = 0; j < ny; ++j)
		for (size_t i = 0; i < nx + 1; ++i)
			u[i + j * (nx + 1)] = 0;
}

void init_v() {
	v.resize(nx * (ny + 1));
	// fill v
	for (size_t j = 0; j < ny + 1; ++j)
		for (size_t i = 0; i < nx; ++i)
			v[i + j * nx] = 0;
}

void init_p() {
	p.resize(nx * ny);
	// fill p
	for (size_t j = 0; j < ny; ++j)
		for (size_t i = 0; i < nx; ++i)
			p[i + j * nx] = 0.0;

	//p[0] = 0;
}

void init_T() {
	T.resize((nx + 1) * (ny + 1));
	// fill T
	for (size_t j = 0; j < ny + 1; ++j)
		for (size_t i = 0; i < nx + 1; ++i)
			T[i + j * (nx + 1)] = 0.5;

	//T[0] = 0.5;
}

void init_C() {
	C.resize((nx + 1) * (ny + 1));
	// fill C
	for (size_t j = 0; j < ny + 1; ++j)
		for (size_t i = 0; i < nx + 1; ++i)
			C[i + j * (nx + 1)] = 0.5;

}

void create_mesh() {
	read_data();
	size_t ix = nx + 1, iy = ny + 1;

	nodes_x.resize(ix);
	step_x = (right_x - left_x) / nx;
	for (size_t i = 0; i < ix; ++i)
		nodes_x[i] = left_x + i * step_x;

	h_x.resize(nx);
	for (size_t i = 0; i < nx; ++i)
		h_x[i] = nodes_x[i + 1] - nodes_x[i];

	bar_x.resize(ix);
	bar_x[0] = h_x[0] / 2;
	bar_x[nx] = h_x[nx - 1] / 2;
	for (size_t i = 1; i < nx; ++i)
		bar_x[i] = (h_x[i] + h_x[i - 1]) / 2;

	nodes_y.resize(iy);
	step_y = (right_y - left_y) / ny;
	for (size_t i = 0; i < iy; ++i)
		nodes_y[i] = left_y + i * step_y;

	h_y.resize(ny);
	for (size_t i = 0; i < ny; ++i)
		h_y[i] = nodes_y[i + 1] - nodes_y[i];

	bar_y.resize(iy);
	bar_y[0] = h_y[0] / 2;
	bar_y[ny] = h_y[ny - 1] / 2;
	for (size_t i = 1; i < ny; ++i)
		bar_y[i] = (h_y[i] + h_y[i - 1]) / 2;

	std::cout << "\nnumber of h_x intervals:\t" << h_x.size();
	std::cout << "\nnumber of h_y intervals:\t" << h_y.size();
	std::cout << "\nnumber of b_x intervals:\t" << bar_x.size();
	std::cout << "\nnumber of b_y intervals:\t" << bar_y.size();
	std::cout << "\nnumber of n_x intervals:\t" << nodes_x.size();
	std::cout << "\nnumber of n_y intervals:\t" << nodes_y.size();

	init_u();
	init_v();
	init_p();
	init_T();
	init_C();
}


//---------------------------------------------------------------------------------
// КОЭФФИЦИЕНТЫ ДЛЯ УРАВНЕНИЯ du/dt

// u_i_(j-1/2)
double A1_u(size_t i, size_t j) {

	// regular nodes
	double res = 0;
	res = -tau * nu * bar_x[i] / bar_y[j];

	return res;
}

// u_(i+1)_(j+1/2)
double A5_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * nu * h_y[j] / h_x[i];
	return res;
}

// u_i_(j+1/2)
double A6_u(size_t i, size_t j)
{
	double res = 0;
	res += bar_x[i] * h_y[j];
	res += tau * nu * h_y[j] / h_x[i];
	res += tau * nu * h_y[j] / h_x[i - 1];
	res += tau * nu * bar_x[i] / bar_y[j + 1];
	res += tau * nu * bar_x[i] / bar_y[j];
	return res;
}

// u_i_(j+1/2)
// учитывая, что du/dy = Ma DT/dx
double A6_u_new(size_t i, size_t j)
{
	double res = 0;
	res += bar_x[i] * h_y[j];
	res += tau * nu * h_y[j] / h_x[i];
	res += tau * nu * h_y[j] / h_x[i - 1];
	//res += tau * nu * bar_x[i] / bar_y[j + 1];
	res += tau * nu * bar_x[i] / bar_y[j];
	return res;
}

// u_(i-1)_(j+1/2)
double A7_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * nu * h_y[j] / h_x[i - 1];
	return res;
}

// u_i_(j+3/2)
double A9_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * nu * bar_x[i] / bar_y[j + 1];
	return res;
}

// v_(i+1/2)_j
double B4_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * h_y[j] * h_x[i] / 4.0;
	// w_i_j :
	//res *= (v[i + j * nx] + v[i - 1 + j * nx]) / bar_x[i] - (u[i + j * (nx + 1)] + u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + j * nx] - v[i - 1 + j * nx]) / bar_x[i] - (u[i + j * (nx + 1)] - u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	return res;
}

// v_(i-1/2)_j
double B5_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * h_y[j] * h_x[i - 1] / 4.0;
	//res = h_y[j] * h_x[i - 1] / 4;
	// w_i_j :
	//res *= (v[i + j * nx] + v[i - 1 + j * nx]) / bar_x[i]- (u[i + j * (nx + 1)] + u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + j * nx] - v[i - 1 + j * nx]) / bar_x[i] - (u[i + j * (nx + 1)] - u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	return res;
}

// v_(i-1/2)_(j+1)
double B6_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * h_y[j] * h_x[i - 1] / 4.0;
	//res = h_y[j] * h_x[i - 1] / 4;
	// w_i_(j+1) :
	//res *= (v[i + (j + 1) * nx] + v[i - 1 + (j + 1) * nx]) / bar_x[i] - (u[i + (j + 1) * (nx + 1)] + u[i + j * (nx + 1)]) / bar_y[j + 1];
	res *= (v[i + (j + 1) * nx] - v[i - 1 + (j + 1) * nx]) / bar_x[i] - (u[i + (j + 1) * (nx + 1)] - u[i + j * (nx + 1)]) / bar_y[j + 1];
	return res;
}

// v_(i+1/2)_(j+1)
double B7_u(size_t i, size_t j)
{
	double res = 0;
	res = -tau * h_y[j] * h_x[i] / 4.0;
	//res = h_y[j] * h_x[i] / 4;
	// w_i_(j+1) :
	//res *= (v[i + (j + 1) * nx] + v[i - 1 + (j + 1) * nx]) / bar_x[i]- (u[i + (j + 1) * (nx + 1)] + u[i + j * (nx + 1)]) / bar_y[j + 1];
	res *= (v[i + (j + 1) * nx] - v[i - 1 + (j + 1) * nx]) / bar_x[i] - (u[i + (j + 1) * (nx + 1)] - u[i + j * (nx + 1)]) / bar_y[j + 1];
	return res;
}

// p_(i+1/2)_(j+1/2)
double C4_u(size_t i, size_t j)
{
	double res = tau * h_y[j];
	return res;
}

// p_(i-1/2)_(j+1/2)
double C5_u(size_t i, size_t j)
{
	double res = -tau * h_y[j];
	return res;
}

//---------------------------------------------------------------------------------
// КОЭФФИЦИЕНТЫ ДЛЯ УРАВНЕНИЯ dv/dt

// v_(i+1/2)_(j-1)
double E1_v(size_t i, size_t j)
{
	double res = 0;
	res = -tau * nu * h_x[i] / h_y[j - 1];
	return res;
}

// v_(i-1/2)_j
double E5_v(size_t i, size_t j)
{
	double res = -tau * nu * bar_y[j] / bar_x[i];
	return res;
}

// v_(i+1/2)_j
double E4_v(size_t i, size_t j)
{
	double res = 0;
	res = bar_y[j] * h_x[i];
	res += tau * nu * bar_y[j] / bar_x[i + 1];
	res += tau * nu * bar_y[j] / bar_x[i];
	res += tau * nu * h_x[i] / h_y[j];
	res += tau * nu * h_x[i] / h_y[j - 1];
	return res;
}

// v_(i+3/2)_j
double E3_v(size_t i, size_t j)
{
	double res = -tau * nu * bar_y[j] / bar_x[i + 1];
	return res;
}

// v_(i+1/2)_(j+1)
double E7_v(size_t i, size_t j)
{
	double res = -tau * nu * h_x[i] / h_y[j];
	return res;
}

// u_i_(j-1/2)
double D1_v(size_t i, size_t j)
{
	double res = tau * h_x[i] * h_y[j - 1] / 4.0;
	// w_i_j :
	//res *= (v[i + j * nx] + v[i - 1 + j * nx]) / bar_x[i]- (u[i + j * (nx + 1)] + u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + j * nx] - v[i - 1 + j * nx]) / bar_x[i] - (u[i + j * (nx + 1)] - u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	return res;
}

// u_(i+1)_(j-1/2)
double D2_v(size_t i, size_t j)
{
	double res = tau * h_x[i] * h_y[j - 1] / 4.0;
	// w_(i+1)_j :
	//res *= (v[i + 1 + j * nx] + v[i + j * nx]) / bar_x[i + 1]	- (u[i + 1 + j * (nx + 1)] + u[i + 1 + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + 1 + j * nx] - v[i + j * nx]) / bar_x[i + 1] - (u[i + 1 + j * (nx + 1)] - u[i + 1 + (j - 1) * (nx + 1)]) / bar_y[j];

	return res;
}

// u_i_(j+1/2)
double D6_v(size_t i, size_t j)
{
	double res = tau * h_y[j] * h_x[i] / 4.0;
	// w_i_j :
	//res *= (v[i + j * nx] + v[i - 1 + j * nx]) / bar_x[i]		- (u[i + j * (nx + 1)] + u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + j * nx] - v[i - 1 + j * nx]) / bar_x[i] - (u[i + j * (nx + 1)] - u[i + (j - 1) * (nx + 1)]) / bar_y[j];
	return res;
}

// u_(i+1)_(j+1/2)
double D5_v(size_t i, size_t j)
{
	double res = tau * h_x[i] * h_y[j] / 4.0;
	// w_(i+1)_j :
	//res *= (v[i + 1 + j * nx] + v[i + j * nx]) / bar_x[i + 1] - (u[i + 1 + j * (nx + 1)] + u[i + 1 + (j - 1) * (nx + 1)]) / bar_y[j];
	res *= (v[i + 1 + j * nx] - v[i + j * nx]) / bar_x[i + 1] - (u[i + 1 + j * (nx + 1)] - u[i + 1 + (j - 1) * (nx + 1)]) / bar_y[j];
	return res;
}

// p_(i+1/2)_(j-1/2)
double F1_v(size_t i, size_t j)
{
	double res = -tau * h_x[i];
	return res;
}

// p_(i+1/2)_(j+1/2)
double F4_v(size_t i, size_t j)
{
	double res = tau * h_x[i];
	return res;
}

//---------------------------------------------------------------------------------
// КОЭФФИЦИЕНТЫ ДЛЯ УРАВНЕНИЯ du/dx + dv/dy = 0

// u_(i+1)_(j+1/2)
double coef0_BLOCK_P(size_t i, size_t j)
{
	//double res = 1 / h_x[i];
	double res = h_y[j];
	return res;
}

// u_i_(j+1/2)
double coef1_BLOCK_P(size_t i, size_t j)
{
	//double res = -1 / h_x[i];
	double res = -h_y[j];
	return res;
}

// v_(i+1/2)_(j+1)
double coef2_BLOCK_P(size_t i, size_t j)
{
	//double res = 1 / h_y[j];
	double res = h_x[i];
	return res;
}

// v_(i+1/2)_j
double coef3_BLOCK_P(size_t i, size_t j)
{
	//double res = -1 / h_y[j];
	double res = -h_x[i];
	return res;
}

// ---------------------------------------------------------------------------------

// matrix
void matrix() {
	// u block
	/*
		matrix_size (nx + 1) * ny + nx * (ny + 1) + nx * ny
					|___________|   |___________|   |_____|
						len_u           len_v        len_p
	*/

	size_t counter_elems = 0;
	size_t counter_strng = 0;

	size_t len_u = (nx + 1) * ny;
	size_t len_v = nx * (ny + 1);
	size_t len_p = nx * ny;
	size_t len_m = len_u + len_v + len_p;

	// j = 0
	for (size_t i = 0; i < nx + 1; ++i) {
		if (i == 0 || i == nx) {
			c_r[counter_strng] = counter_elems + 1;
			c_c[counter_elems] = i + 1;
			c_m[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else { // check filling columns
			c_r[counter_strng] = counter_elems + 1;
			c_c[counter_elems] = i - 1 + 1;
			c_m[counter_elems] = A7_u(i, 0); // u_0_1/2
			++counter_elems;
			c_c[counter_elems] = i + 1;
			c_m[counter_elems] = A6_u(i, 0); // u_1_1/2
			++counter_elems;
			c_c[counter_elems] = i + 1 + 1;
			c_m[counter_elems] = A5_u(i, 0); // u_2_1/2
			++counter_elems;
			c_c[counter_elems] = i + nx + 1 + 1;
			c_m[counter_elems] = A9_u(i, 0); // u_1_3/2
			++counter_elems;

			c_c[counter_elems] = len_u + nx + i - 1 + 1;
			c_m[counter_elems] = B6_u(i, 0); // v_(i-1/2)_(j+1)
			++counter_elems;
			c_c[counter_elems] = len_u + nx + i + 1;
			c_m[counter_elems] = B7_u(i, 0); // v_(i+1/2)_(j+1)
			++counter_elems;

			c_c[counter_elems] = i - 1 + len_u + len_v + 1;
			c_m[counter_elems] = C5_u(i, 0); // p_(i-1/2)_(j+1/2)
			++counter_elems;
			c_c[counter_elems] = i + len_u + len_v + 1;
			c_m[counter_elems] = C4_u(i, 0); // p_(i+1/2)_(j+1/2)
			++counter_elems;
			++counter_strng;
		}
	}

	// j from 1 to ny - 2
	for (size_t j = 1; j < ny - 1; ++j) {
		for (size_t i = 0; i < nx + 1; ++i) {
			if (i == 0 || i == nx) {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = 1.0;
				++counter_elems;
				++counter_strng;
			}
			else {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = A1_u(i, j); // u_i_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i - 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = A7_u(i, j); // u_(i-1)_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = A6_u(i, j); // u_i_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = A5_u(i, j); // u_(i+1)_(j+1/2) 
				++counter_elems;
				c_c[counter_elems] = i + (j + 1) * (nx + 1) + 1;
				c_m[counter_elems] = A9_u(i, j); // u_i_(j+3/2)
				++counter_elems;

				c_c[counter_elems] = i - 1 + len_u + j * nx + 1;
				c_m[counter_elems] = B5_u(i, j); // v_(i-1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + len_u + j * nx + 1;
				c_m[counter_elems] = B4_u(i, j); // v_(i+1/2)_j
				++counter_elems;
				c_c[counter_elems] = len_u + i - 1 + (j + 1) * nx + 1;
				//c_m[counter_elems] = B6_u(i, 0); // v_(i-1/2)_(j+1)
				c_m[counter_elems] = B6_u(i, j); // v_(i-1/2)_(j+1)
				++counter_elems;
				c_c[counter_elems] = len_u + i + (j + 1) * nx + 1;
				//c_m[counter_elems] = B7_u(i, 0); // v_(i+1/2)_(j+1)
				c_m[counter_elems] = B7_u(i, j); // v_(i+1/2)_(j+1)
				++counter_elems;

				c_c[counter_elems] = i - 1 + j * nx + len_u + len_v + 1;
				c_m[counter_elems] = C5_u(i, j); // p_(i-1/2)_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + len_v + 1;
				c_m[counter_elems] = C4_u(i, j); // p_(i+1/2)_(j+1/2)
				++counter_elems;
				++counter_strng;
			}
		}
	}

	// j = ny - 1
	for (size_t i = 0; i < nx + 1; ++i) {
		if (i == 0 || i == nx) {
			c_r[counter_strng] = counter_elems + 1;
			c_c[counter_elems] = i + (ny - 1) * (nx + 1) + 1;
			c_m[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else { // check filling columns
			c_r[counter_strng] = counter_elems + 1;
			c_c[counter_elems] = i + (ny - 2) * (nx + 1) + 1;
			c_m[counter_elems] = A1_u(i, ny - 1); // u_i_(j-1/2)
			++counter_elems;
			c_c[counter_elems] = i - 1 + (ny - 1) * (nx + 1) + 1;
			c_m[counter_elems] = A7_u(i, ny - 1); // u_(i-1)_(j+1/2)
			++counter_elems;

			c_c[counter_elems] = i + (ny - 1) * (nx + 1) + 1;
			//c_m[counter_elems] = A6_u(i, (ny - 1)); // u_i_(j+1/2)
			// учитывая, что du/dy = Ma dT/dx:
			c_m[counter_elems] = A6_u_new(i, (ny - 1)); // u_i_(j+1/2)
			++counter_elems;
			c_c[counter_elems] = i + 1 + (ny - 1) * (nx + 1) + 1;
			c_m[counter_elems] = A5_u(i, (ny - 1)); // u_(i+1)_(j+1/2) 
			++counter_elems;

			c_c[counter_elems] = i - 1 + len_u + (ny - 1) * nx + 1;
			c_m[counter_elems] = B5_u(i, (ny - 1)); // v_(i-1/2)_j
			++counter_elems;
			c_c[counter_elems] = i + len_u + (ny - 1) * nx + 1;
			c_m[counter_elems] = B4_u(i, (ny - 1)); // v_(i+1/2)_j
			++counter_elems;

			c_c[counter_elems] = i - 1 + (ny - 1) * nx + len_u + len_v + 1;
			c_m[counter_elems] = C5_u(i, (ny - 1)); // p_(i-1/2)_(j+1/2)
			++counter_elems;
			c_c[counter_elems] = i + (ny - 1) * nx + len_u + len_v + 1;
			c_m[counter_elems] = C4_u(i, (ny - 1)); // p_(i+1/2)_(j+1/2)
			++counter_elems;
			++counter_strng;
		}
	}

	// v block
	// j = 0

	for (size_t i = 0; i < nx; ++i) {
		c_r[counter_strng] = counter_elems + 1;
		c_c[counter_elems] = i + len_u + 1;
		c_m[counter_elems] = 1.0;
		++counter_elems;
		++counter_strng;
	}

	// j from 1 to ny - 1
	for (size_t j = 1; j < ny; ++j) {
		for (size_t i = 0; i < nx; ++i) {
			if (i == 0) {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D1_v(i, j) * 0; // u_i_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D2_v(i, j); // u_(i+1)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = D6_v(i, j) * 0; // u_i_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = D5_v(i, j); // u_(i+1)_(j+1/2)
				++counter_elems;


				c_c[counter_elems] = i + (j - 1) * nx + len_u + 1;
				c_m[counter_elems] = E1_v(i, j); // v_(i+1/2)_(j-1)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + 1;
				c_m[counter_elems] = E4_v(i, j); // v_(i+1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * nx + len_u + 1;
				c_m[counter_elems] = E3_v(i, j); // v_(i+3/2)_j
				++counter_elems;
				c_c[counter_elems] = i + (j + 1) * nx + len_u + 1;
				c_m[counter_elems] = E7_v(i, j); // v_(i+1/2)_(j+1)
				++counter_elems;

				c_c[counter_elems] = i + (j - 1) * nx + len_u + len_v + 1;
				c_m[counter_elems] = F1_v(i, j); // p_(i+1/2)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + len_v + 1;
				c_m[counter_elems] = F4_v(i, j); // p_(i+1/2)_(j+1/2)
				++counter_elems;
				++counter_strng;
			}
			else if (i == nx - 1) {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D1_v(i, j); // u_i_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D2_v(i, j) * 0; // u_(i+1)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = D6_v(i, j); // u_i_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = D5_v(i, j) * 0; // u_(i+1)_(j+1/2)
				++counter_elems;


				c_c[counter_elems] = i + (j - 1) * nx + len_u + 1;
				c_m[counter_elems] = E1_v(i, j); // v_(i+1/2)_(j-1)
				++counter_elems;
				c_c[counter_elems] = i - 1 + j * nx + len_u + 1;
				c_m[counter_elems] = E5_v(i, j); // v_(i-1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + 1;
				c_m[counter_elems] = E4_v(i, j); // v_(i+1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + (j + 1) * nx + len_u + 1;
				c_m[counter_elems] = E7_v(i, j); // v_(i+1/2)_(j+1)
				++counter_elems;

				c_c[counter_elems] = i + (j - 1) * nx + len_u + len_v + 1;
				c_m[counter_elems] = F1_v(i, j); // p_(i+1/2)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + len_v + 1;
				c_m[counter_elems] = F4_v(i, j); // p_(i+1/2)_(j+1/2)
				++counter_elems;
				++counter_strng;
			}
			else {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D1_v(i, j); // u_i_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + (j - 1) * (nx + 1) + 1;
				c_m[counter_elems] = D2_v(i, j); // u_(i+1/2)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = D6_v(i, j); // u_i_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = D5_v(i, j); // u_(i+1)_(j+1/2)
				++counter_elems;


				c_c[counter_elems] = i + (j - 1) * nx + len_u + 1;
				c_m[counter_elems] = E1_v(i, j); // v_(i+1/2)_(j-1)
				++counter_elems;
				c_c[counter_elems] = i - 1 + j * nx + len_u + 1;
				c_m[counter_elems] = E5_v(i, j); // v_(i-1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + 1;
				c_m[counter_elems] = E4_v(i, j); // v_(i+1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * nx + len_u + 1;
				c_m[counter_elems] = E3_v(i, j); // v_(i+3/2)_j
				++counter_elems;
				c_c[counter_elems] = i + (j + 1) * nx + len_u + 1;
				c_m[counter_elems] = E7_v(i, j); // v_(i+1/2)_(j+1)
				++counter_elems;

				c_c[counter_elems] = i + (j - 1) * nx + len_u + len_v + 1;
				c_m[counter_elems] = F1_v(i, j); // p_(i+1/2)_(j-1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + len_v + 1;
				c_m[counter_elems] = F4_v(i, j); // p_(i+1/2)_(j+1/2)
				++counter_elems;
				++counter_strng;
			}
		}
	}

	// j = ny
	for (size_t i = 0; i < nx; ++i) {
		c_r[counter_strng] = counter_elems + 1;
		c_c[counter_elems] = i + ny * nx + len_u + 1;
		c_m[counter_elems] = 1.0;
		++counter_elems;
		++counter_strng;
	}

	// p block

	for (size_t j = 0; j < ny; ++j) {
		for (size_t i = 0; i < nx; ++i) {
			if ((i == 0) && (j == 0)) {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = len_u + len_v + 1;
				c_m[counter_elems] = 1.0;
				++counter_elems;
				++counter_strng;
			}
			else {
				c_r[counter_strng] = counter_elems + 1;
				c_c[counter_elems] = i + j * (nx + 1) + 1;
				c_m[counter_elems] = coef1_BLOCK_P(i, j); // u_i_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m[counter_elems] = coef0_BLOCK_P(i, j); // u_(i+1)_(j+1/2)
				++counter_elems;
				c_c[counter_elems] = i + j * nx + len_u + 1;
				c_m[counter_elems] = coef3_BLOCK_P(i, j); // v_(i+1/2)_j
				++counter_elems;
				c_c[counter_elems] = i + (j + 1) * nx + len_u + 1;
				c_m[counter_elems] = coef2_BLOCK_P(i, j); // v_(i+1/2)_(j+1)
				++counter_elems;
				++counter_strng;
			}
		}

	} // p block 

	c_r[counter_strng] = counter_elems + 1;


	/*std::cout << "\n\nc_r\t size:" << size_c+1<<" "<<counter_strng+1 << " elems\n";
	for (size_t i = 0; i <= size_c; ++i)
		std::cout << c_r[i] << " ";

			std::cout << "\n\nc_c \t size:" << counter_elems << " elems\n";
		for (size_t i = 0; i < size_m; ++i)
			std::cout << c_c[i] << " ";



		std::cout << "\n\nc_m\t size:" << counter_elems << " elems\n";
		for (size_t i = 0; i < size_m; ++i)
			std::cout << c_m[i] << " ";*/


} // matrix


void matrix_T()
{
	/*
		matrix_size (nx + 1) * (ny + 1)
	*/

	size_t counter_elems = 0;
	size_t counter_strng = 0;

	// j = 0
	for (size_t i = 0; i < nx + 1; ++i) {
		if ((i == 0) || (i == nx)) {
			c_r_T[counter_strng] = counter_elems + 1;
			c_c_T[counter_elems] = i + 1;
			c_m_T[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else
		{ // check filling columns

			// T_i-1_j
			c_r_T[counter_strng] = counter_elems + 1;
			c_c_T[counter_elems] = i - 1 + 1;
			c_m_T[counter_elems] = (-
				0.125 * h_y[0] * u[(i - 1)] -
				0.125 * h_y[0] * u[i] -
				bar_y[0] / (pr * h_x[i - 1]))
				* tau;
			++counter_elems;

			// T_i_j
			c_c_T[counter_elems] = i + 1;

			c_m_T[counter_elems] = bar_y[0] * bar_x[i] +
				tau * (
					bar_y[0] / (pr * h_x[i]) +
					bar_y[0] / (pr * h_x[i - 1]) +
					bar_x[i] / (pr * h_y[0]) +
					0.125 * h_y[0] * u[i + 1] -
					0.125 * h_y[0] * u[i - 1] +
					0.125 * v[i + 1 * nx] * h_x[i] +
					0.125 * v[i - 1 + 1 * nx] * h_x[i - 1]
					);

			++counter_elems;

			// T_i+1_j
			c_c_T[counter_elems] = i + 1 + 1;
			c_m_T[counter_elems] = tau * (
				0.125 * h_y[0] * u[(i + 1) + 0 * (nx + 1)] +
				0.125 * h_y[0] * u[i + 0 * (nx + 1)] -
				bar_y[0] / (pr * h_x[i])
				);
			++counter_elems;

			// T_i_j+1
			c_c_T[counter_elems] = i + (nx + 1) + 1;
			c_m_T[counter_elems] = tau * (
				0.125 * v[i + 0 * nx] * h_x[i] +
				0.125 * v[i - 1 + 0 * nx] * h_x[i - 1] +
				0.125 * v[i + 1 + 1 * nx] * h_x[i] +
				0.125 * v[i - 1 + 1 * nx] * h_x[i - 1] -
				bar_x[i] / (pr * h_y[0])
				);
			++counter_elems;

			++counter_strng;
		}
	}

	// j from 1 to ny - 1
	for (size_t j = 1; j < ny; ++j) {
		for (size_t i = 0; i < nx + 1; ++i) {
			if ((i == 0) || (i == nx)) {
				c_r_T[counter_strng] = counter_elems + 1;
				c_c_T[counter_elems] = i + j * (nx + 1) + 1;
				c_m_T[counter_elems] = 1.0;
				++counter_elems;
				++counter_strng;
			}
			else
			{
				c_r_T[counter_strng] = counter_elems + 1;

				//T_i_j-1
				c_c_T[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m_T[counter_elems] = tau * (-
					0.125 * v[i + (j - 1) * nx] * h_x[i] -
					0.125 * v[i - 1 + (j - 1) * nx] * h_x[i - 1] -
					0.125 * v[i + j * nx] * h_x[i] -
					0.125 * v[i - 1 + j * nx] * h_x[i - 1] -
					bar_x[i] / (pr * h_y[j - 1]));
				++counter_elems;

				// T_i-1_j
				c_c_T[counter_elems] = i - 1 + j * (nx + 1) + 1;
				c_m_T[counter_elems] = (-
					0.125 * h_y[j] * u[(i - 1) + j * (nx + 1)] -
					0.125 * h_y[j - 1] * u[(i - 1) + (j - 1) * (nx + 1)] -
					0.125 * h_y[j] * u[i + j * (nx + 1)] -
					0.125 * h_y[j - 1] * u[i + (j - 1) * (nx + 1)] -
					bar_y[j] / (pr * h_x[i - 1])) * tau;
				++counter_elems;

				// T_i_j
				c_c_T[counter_elems] = i + j * (nx + 1) + 1;
				c_m_T[counter_elems] = bar_y[j] * bar_x[i] +
					tau * (
						bar_y[j] / (pr * h_x[i]) +
						bar_y[j] / (pr * h_x[i - 1]) +
						bar_x[i] / (pr * h_y[j]) +
						bar_x[i] / (pr * h_y[j - 1]) +
						0.125 * h_y[j] * u[i + 1 + j * (nx + 1)] +
						0.125 * h_y[j - 1] * u[i + 1 + (j - 1) * (nx + 1)] -
						0.125 * h_y[j] * u[i - 1 + j * (nx + 1)] -
						0.125 * h_y[j - 1] * u[i - 1 + (j - 1) * (nx + 1)] +
						0.125 * h_x[i] * v[i + (j + 1) * nx] +
						0.125 * h_x[i - 1] * v[i - 1 + (j + 1) * nx] -
						0.125 * h_x[i] * v[i + (j - 1) * nx] -
						0.125 * h_x[i - 1] * v[i - 1 + (j - 1) * nx]
						);
				++counter_elems;

				// T_i+1_j
				c_c_T[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m_T[counter_elems] = tau * (
					0.125 * h_y[j] * u[i + 1 + j * (nx + 1)] +
					0.125 * h_y[j - 1] * u[i + 1 + (j - 1) * (nx + 1)] +
					0.125 * h_y[j] * u[i + j * (nx + 1)] +
					0.125 * h_y[j - 1] * u[i + (j - 1) * (nx + 1)] -
					bar_y[j] / (pr * h_x[i]));
				++counter_elems;

				// T_i_j+1
				c_c_T[counter_elems] = i + (j + 1) * (nx + 1) + 1;
				c_m_T[counter_elems] = tau * (
					0.125 * v[i + j * nx] * h_x[i] +
					0.125 * v[i - 1 + j * nx] * h_x[i - 1] +
					0.125 * v[i + (j + 1) * nx] * h_x[i] +
					0.125 * v[i - 1 + (j + 1) * nx] * h_x[i - 1] -
					bar_x[i] / (pr * h_y[j]));
				++counter_elems;

				++counter_strng;
			}
		}
	}

	// j = ny
	for (size_t i = 0; i < nx + 1; ++i) {
		if ((i == 0) || (i == nx)) {
			c_r_T[counter_strng] = counter_elems + 1;
			c_c_T[counter_elems] = i + ny * (nx + 1) + 1;
			c_m_T[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else
		{ // check filling columns
			c_r_T[counter_strng] = counter_elems + 1;

			//T_i_j-1
			c_c_T[counter_elems] = i + (ny - 1) * (nx + 1) + 1;
			c_m_T[counter_elems] = tau * (-
				0.125 * v[i + (ny - 1) * nx] * h_x[i] -
				0.125 * v[i - 1 + (ny - 1) * nx] * h_x[i - 1] -
				0.125 * v[i + ny * nx] * h_x[i] -
				0.125 * v[i - 1 + ny * nx] * h_x[i - 1]
				- bar_x[i] / (pr * h_y[ny - 1]));
			++counter_elems;

			// T_i-1_j
			c_c_T[counter_elems] = i - 1 + ny * (nx + 1) + 1;
			c_m_T[counter_elems] = (-
				0.125 * h_y[ny - 1] * u[(i - 1) + (ny - 1) * (nx + 1)] -
				0.125 * h_y[ny - 1] * u[i + (ny - 1) * (nx + 1)] -
				bar_y[ny] / (pr * h_x[i - 1])) * tau;
			++counter_elems;

			// T_i_j
			c_c_T[counter_elems] = i + ny * (nx + 1) + 1;
			c_m_T[counter_elems] = bar_y[ny] * bar_x[i] +
				tau * (
					bar_y[ny] / (pr * h_x[i]) +
					bar_y[ny] / (pr * h_x[i - 1]) +
					bar_x[i] / (pr * h_y[ny - 1]) +
					0.125 * h_y[ny - 1] * u[i + 1 + (ny - 1) * (nx + 1)] -
					0.125 * h_y[ny - 1] * u[i - 1 + (ny - 1) * (nx + 1)] -
					0.125 * h_x[i] * v[i + (ny - 1) * nx] -
					0.125 * h_x[i - 1] * v[i - 1 + (ny - 1) * nx]
					);
			++counter_elems;

			// T_i+1_j
			c_c_T[counter_elems] = i + 1 + ny * (nx + 1) + 1;
			c_m_T[counter_elems] = tau * (
				0.125 * h_y[ny - 1] * u[i + 1 + (ny - 1) * (nx + 1)] +
				0.125 * h_y[ny - 1] * u[i + (ny - 1) * (nx + 1)] -
				bar_y[ny] / (pr * h_x[i]));
			++counter_elems;

			++counter_strng;
		}
	}

	c_r_T[counter_strng] = counter_elems + 1;


	/*std::cout << "\nc_c_T: \n";
	for (int i = 0; i < size_m_T;i++)
		std::cout << c_c_T[i] << " ";
	std::cout << "\n\n";

	std::cout << "\nc_m_T: \n";
	for (int i = 0; i < size_m_T;i++)
		std::cout << c_m_T[i] << " ";
	std::cout << "\n\n";

	std::cout << "\nc_r_T: \n";
	for (int i = 0; i <= size_c_T;i++)
		std::cout << c_r_T[i] << " ";
	std::cout << "\n\n";*/


}

void matrix_C()
{
	/*
		matrix_size (nx + 1) * (ny + 1)
	*/

	size_t counter_elems = 0;
	size_t counter_strng = 0;

	// j = 0
	for (size_t i = 0; i < nx + 1; ++i) {
		if ((i == 0) || (i == nx)) {
			c_r_C[counter_strng] = counter_elems + 1;
			c_c_C[counter_elems] = i + 1;
			c_m_C[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else
		{ // check filling columns

			// C_i-1_j
			c_r_C[counter_strng] = counter_elems + 1;
			c_c_C[counter_elems] = i - 1 + 1;
			c_m_C[counter_elems] = (-
				0.125 * h_y[0] * u[(i - 1)] -
				0.125 * h_y[0] * u[i] -
				bar_y[0] / (Sc * h_x[i - 1]))
				* tau;
			++counter_elems;

			// C_i_j
			c_c_C[counter_elems] = i + 1;

			c_m_C[counter_elems] = bar_y[0] * bar_x[i] +
				tau * (
					bar_y[0] / (Sc * h_x[i]) +
					bar_y[0] / (Sc * h_x[i - 1]) +
					bar_x[i] / (Sc * h_y[0]) +
					0.125 * h_y[0] * u[i + 1] -
					0.125 * h_y[0] * u[i - 1] +
					0.125 * v[i + 1 * nx] * h_x[i] +
					0.125 * v[i - 1 + 1 * nx] * h_x[i - 1]
					);

			++counter_elems;

			// C_i+1_j
			c_c_C[counter_elems] = i + 1 + 1;
			c_m_C[counter_elems] = tau * (
				0.125 * h_y[0] * u[(i + 1) + 0 * (nx + 1)] +
				0.125 * h_y[0] * u[i + 0 * (nx + 1)] -
				bar_y[0] / (Sc * h_x[i])
				);
			++counter_elems;

			// C_i_j+1
			c_c_C[counter_elems] = i + (nx + 1) + 1;
			c_m_C[counter_elems] = tau * (
				0.125 * v[i + 0 * nx] * h_x[i] +
				0.125 * v[i - 1 + 0 * nx] * h_x[i - 1] +
				0.125 * v[i + 1 + 1 * nx] * h_x[i] +
				0.125 * v[i - 1 + 1 * nx] * h_x[i - 1] -
				bar_x[i] / (Sc * h_y[0])
				);
			++counter_elems;

			++counter_strng;
		}
	}

	// j from 1 to ny - 1
	for (size_t j = 1; j < ny; ++j) {
		for (size_t i = 0; i < nx + 1; ++i) {
			if ((i == 0) || (i == nx)) {
				c_r_C[counter_strng] = counter_elems + 1;
				c_c_C[counter_elems] = i + j * (nx + 1) + 1;
				c_m_C[counter_elems] = 1.0;
				++counter_elems;
				++counter_strng;
			}
			else
			{
				c_r_C[counter_strng] = counter_elems + 1;

				//C_i_j-1
				c_c_C[counter_elems] = i + (j - 1) * (nx + 1) + 1;
				c_m_C[counter_elems] = tau * (-
					0.125 * v[i + (j - 1) * nx] * h_x[i] -
					0.125 * v[i - 1 + (j - 1) * nx] * h_x[i - 1] -
					0.125 * v[i + j * nx] * h_x[i] -
					0.125 * v[i - 1 + j * nx] * h_x[i - 1] -
					bar_x[i] / (Sc * h_y[j - 1]));
				++counter_elems;

				// C_i-1_j
				c_c_C[counter_elems] = i - 1 + j * (nx + 1) + 1;
				c_m_C[counter_elems] = (-
					0.125 * h_y[j] * u[(i - 1) + j * (nx + 1)] -
					0.125 * h_y[j - 1] * u[(i - 1) + (j - 1) * (nx + 1)] -
					0.125 * h_y[j] * u[i + j * (nx + 1)] -
					0.125 * h_y[j - 1] * u[i + (j - 1) * (nx + 1)] -
					bar_y[j] / (Sc * h_x[i - 1])) * tau;
				++counter_elems;

				// C_i_j
				c_c_C[counter_elems] = i + j * (nx + 1) + 1;
				c_m_C[counter_elems] = bar_y[j] * bar_x[i] +
					tau * (
						bar_y[j] / (Sc * h_x[i]) +
						bar_y[j] / (Sc * h_x[i - 1]) +
						bar_x[i] / (Sc * h_y[j]) +
						bar_x[i] / (Sc * h_y[j - 1]) +
						0.125 * h_y[j] * u[i + 1 + j * (nx + 1)] +
						0.125 * h_y[j - 1] * u[i + 1 + (j - 1) * (nx + 1)] -
						0.125 * h_y[j] * u[i - 1 + j * (nx + 1)] -
						0.125 * h_y[j - 1] * u[i - 1 + (j - 1) * (nx + 1)] +
						0.125 * h_x[i] * v[i + (j + 1) * nx] +
						0.125 * h_x[i - 1] * v[i - 1 + (j + 1) * nx] -
						0.125 * h_x[i] * v[i + (j - 1) * nx] -
						0.125 * h_x[i - 1] * v[i - 1 + (j - 1) * nx]
						);
				++counter_elems;

				// C_i+1_j
				c_c_C[counter_elems] = i + 1 + j * (nx + 1) + 1;
				c_m_C[counter_elems] = tau * (
					0.125 * h_y[j] * u[i + 1 + j * (nx + 1)] +
					0.125 * h_y[j - 1] * u[i + 1 + (j - 1) * (nx + 1)] +
					0.125 * h_y[j] * u[i + j * (nx + 1)] +
					0.125 * h_y[j - 1] * u[i + (j - 1) * (nx + 1)] -
					bar_y[j] / (Sc * h_x[i]));
				++counter_elems;

				// C_i_j+1
				c_c_C[counter_elems] = i + (j + 1) * (nx + 1) + 1;
				c_m_C[counter_elems] = tau * (
					0.125 * v[i + j * nx] * h_x[i] +
					0.125 * v[i - 1 + j * nx] * h_x[i - 1] +
					0.125 * v[i + (j + 1) * nx] * h_x[i] +
					0.125 * v[i - 1 + (j + 1) * nx] * h_x[i - 1] -
					bar_x[i] / (Sc * h_y[j]));
				++counter_elems;

				++counter_strng;
			}
		}
	}

	// j = ny
	for (size_t i = 0; i < nx + 1; ++i) {
		if ((i == 0) || (i == nx)) {
			c_r_C[counter_strng] = counter_elems + 1;
			c_c_C[counter_elems] = i + ny * (nx + 1) + 1;
			c_m_C[counter_elems] = 1.0;
			++counter_elems;
			++counter_strng;
		}
		else
		{ // check filling columns
			c_r_C[counter_strng] = counter_elems + 1;

			//C_i_j-1
			c_c_C[counter_elems] = i + (ny - 1) * (nx + 1) + 1;
			c_m_C[counter_elems] = tau * (-
				0.125 * v[i + (ny - 1) * nx] * h_x[i] -
				0.125 * v[i - 1 + (ny - 1) * nx] * h_x[i - 1] -
				0.125 * v[i + ny * nx] * h_x[i] -
				0.125 * v[i - 1 + ny * nx] * h_x[i - 1]
				- bar_x[i] / (Sc * h_y[ny - 1]));
			++counter_elems;

			// C_i-1_j
			c_c_C[counter_elems] = i - 1 + ny * (nx + 1) + 1;
			c_m_C[counter_elems] = (-
				0.125 * h_y[ny - 1] * u[(i - 1) + (ny - 1) * (nx + 1)] -
				0.125 * h_y[ny - 1] * u[i + (ny - 1) * (nx + 1)] -
				bar_y[ny] / (Sc * h_x[i - 1])) * tau;
			++counter_elems;

			// C_i_j
			c_c_C[counter_elems] = i + ny * (nx + 1) + 1;
			c_m_C[counter_elems] = bar_y[ny] * bar_x[i] +
				tau * (
					bar_y[ny] / (Sc * h_x[i]) +
					bar_y[ny] / (Sc * h_x[i - 1]) +
					bar_x[i] / (Sc * h_y[ny - 1]) +
					0.125 * h_y[ny - 1] * u[i + 1 + (ny - 1) * (nx + 1)] -
					0.125 * h_y[ny - 1] * u[i - 1 + (ny - 1) * (nx + 1)] -
					0.125 * h_x[i] * v[i + (ny - 1) * nx] -
					0.125 * h_x[i - 1] * v[i - 1 + (ny - 1) * nx]
					);
			++counter_elems;

			// C_i+1_j
			c_c_C[counter_elems] = i + 1 + ny * (nx + 1) + 1;
			c_m_C[counter_elems] = tau * (
				0.125 * h_y[ny - 1] * u[i + 1 + (ny - 1) * (nx + 1)] +
				0.125 * h_y[ny - 1] * u[i + (ny - 1) * (nx + 1)] -
				bar_y[ny] / (Sc * h_x[i]));
			++counter_elems;

			++counter_strng;
		}
	}

	c_r_C[counter_strng] = counter_elems + 1;


	/*std::cout << "\nc_c_T: \n";
	for (int i = 0; i < size_m_T;i++)
		std::cout << c_c_T[i] << " ";
	std::cout << "\n\n";

	std::cout << "\nc_m_T: \n";
	for (int i = 0; i < size_m_T;i++)
		std::cout << c_m_T[i] << " ";
	std::cout << "\n\n";

	std::cout << "\nc_r_T: \n";
	for (int i = 0; i <= size_c_T;i++)
		std::cout << c_r_T[i] << " ";
	std::cout << "\n\n";*/


}


void right_part() {
	size_t counter = 0;
	size_t len_u = (nx + 1) * ny;
	size_t len_v = nx * (ny + 1);
	size_t len_p = nx * ny;

	for (size_t j = 0; j < ny; ++j)
		for (size_t i = 0; i < nx + 1; ++i)
			if (i == 0 || i == nx)
				c_b[counter++] = 0;
			else if (j == ny - 1)
				/*c_b[counter++] = bar_x[i] * h_y[j] * u[i + j * (nx + 1)]
				+ tau * nu * bar_x[i] / bar_y[j + 1] * u_top;*/

				//// учитывая, что du/dy = Ma dT/dx:
				//c_b[counter++] = bar_x[i] * h_y[j] * u[i + j * (nx + 1)]
				//- nu * tau * Ma / pr * (T[i + 1 + ny * (nx + 1)] - T[i - 1 + ny * (nx + 1)]) / 2.0;
	            
	// учитывая, что du/dy = Ma dT/dx + Ma_C dC/dx:
	c_b[counter++] = bar_x[i] * h_y[j] * u[i + j * (nx + 1)]
		- nu * tau * Ma / pr * (T[i + 1 + ny * (nx + 1)] - T[i - 1 + ny * (nx + 1)]) / 2.0 -
				nu * tau * Ma_C / Sc * (C[i + 1 + ny * (nx + 1)] - C[i - 1 + ny * (nx + 1)]) / 2.0;
			else
				c_b[counter++] = bar_x[i] * h_y[j] * u[i + j * (nx + 1)];

	for (size_t j = 0; j < ny + 1; ++j)
		for (size_t i = 0; i < nx; ++i)
			if ((j == 0) || (j == ny))
				c_b[counter++] = 0;
			else
				c_b[counter++] = h_x[i] * bar_y[j] * v[i + j * nx] +
				tau * 0.5 * h_x[i] * bar_y[j] * gr * (T[i + j * (nx + 1)] + T[i + 1 + j * (nx + 1)]) / pr;

	for (size_t j = 0; j < ny; ++j)
		for (size_t i = 0; i < nx; ++i)
			if (i == 0 && j == 0)
				c_b[counter++] = 0.0;
			else
				c_b[counter++] = 0;


	/*std::cout << "\n\n\n\n";
	std::cout << "right_part size: " << size_c << std::endl;
	for (size_t i = 0; i < size_c; ++i)
		std::cout << c_b[i]<<" ";
	std::cout << "\n";*/


}

void right_part_T() {
	size_t counter = 0;

	c_b_T[counter++] = 1.0;
	for (size_t i = 1; i < nx; ++i)
		c_b_T[counter++] = bar_x[i] * bar_y[0] * T[i + 0 * (nx + 1)];
	c_b_T[counter++] = 0;


	for (size_t j = 1; j < ny; ++j)
	{
		c_b_T[counter++] = 1.0;


		for (size_t i = 1; i < nx; ++i)
		{
			c_b_T[counter++] = bar_x[i] * bar_y[j] * T[i + j * (nx + 1)];
		}

		c_b_T[counter++] = 0;



	}

	c_b_T[counter++] = 1.0;
	for (size_t i = 1; i < nx; ++i)
		c_b_T[counter++] = bar_x[i] * bar_y[ny] * T[i + ny * (nx + 1)];
	c_b_T[counter++] = 0;


	/*std::cout << "\n\nright_part_T:\n";
	for (size_t i = 0; i < size_c_T; i++)
		std::cout << c_b_T[i] << " ";
	std::cout << "\n";*/

}

void right_part_C() {
	size_t counter = 0;

	c_b_C[counter++] = 0.0;
	for (size_t i = 1; i < nx; ++i)
		c_b_C[counter++] = bar_x[i] * bar_y[0] * C[i + 0 * (nx + 1)];
	c_b_C[counter++] = 1.0;


	for (size_t j = 1; j < ny; ++j)
	{
		c_b_C[counter++] = 0.0;


		for (size_t i = 1; i < nx; ++i)
		{
			c_b_C[counter++] = bar_x[i] * bar_y[j] * C[i + j * (nx + 1)];
		}

		c_b_C[counter++] = 1.0;



	}

	c_b_C[counter++] = 0.0;
	for (size_t i = 1; i < nx; ++i)
		c_b_C[counter++] = bar_x[i] * bar_y[ny] * C[i + ny * (nx + 1)];
	c_b_C[counter++] = 1.0;


	/*std::cout << "\n\nright_part_T:\n";
	for (size_t i = 0; i < size_c_T; i++)
		std::cout << c_b_T[i] << " ";
	std::cout << "\n";*/

}

void write_file(std::string var, size_t ctr, std::vector<double> vec, std::vector<double> x_step, std::vector<double> y_step, size_t dim_x, size_t dim_y) {
	std::string file_name = PATH + var + "_" + std::to_string(ctr) + ".vtk";
	std::ofstream file;
	file.open(file_name);

	file << "# vtk DataFile Version 2.0\n";
	file << "V.vtk\n";
	file << "ASCII\n";
	file << "DATASET STRUCTURED_GRID\n";
	file << "DIMENSIONS " << dim_x << " " << dim_y << " " << 1 << "\n";
	file << "Points " << dim_x * dim_y << " float\n";

	for (size_t j = 0; j < dim_y; ++j)
		for (size_t i = 0; i < dim_x; ++i)
			file << nodes_x[i] + x_step[i] / 2 << " " << nodes_y[j] + y_step[j] / 2 << " " << 0 << "\n";

	file << "Point_data " << vec.size() << "\n";
	file << "SCALARS " + var + " float 1\n";
	file << "LOOKUP_TABLE default\n";

	for (size_t j = 0; j < dim_y; ++j) {
		for (size_t i = 0; i < dim_x; ++i)
			file << vec[i + j * dim_x] << " ";
		file << "\n\n";
	}
	file.close();
}

void write_file_T(std::string var, size_t ctr, std::vector<double> vec, std::vector<double> x_step, std::vector<double> y_step, size_t dim_x, size_t dim_y) {
	std::string file_name = PATH + var + "_" + std::to_string(ctr) + ".vtk";
	std::ofstream file;
	file.open(file_name);

	file << "# vtk DataFile Version 2.0\n";
	file << "V.vtk\n";
	file << "ASCII\n";
	file << "DATASET STRUCTURED_GRID\n";
	file << "DIMENSIONS " << dim_x << " " << dim_y << " " << 1 << "\n";
	file << "Points " << dim_x * dim_y << " float\n";

	for (size_t j = 0; j < dim_y; ++j)
		for (size_t i = 0; i < dim_x; ++i)
			file << nodes_x[i] << " " << nodes_y[j] << " " << 0 << "\n";

	file << "Point_data " << vec.size() << "\n";
	file << "SCALARS " + var + " float 1\n";
	file << "LOOKUP_TABLE default\n";

	for (size_t j = 0; j < dim_y; ++j) {
		for (size_t i = 0; i < dim_x; ++i)
			file << vec[i + j * dim_x] << " ";
		file << "\n\n";
	}
	file.close();
}

void write_vector_file(std::string var1, std::string var2, size_t ctr, std::vector<double> vec1, std::vector<double> vec2, size_t dim_x, size_t dim_y) {
	std::string file_name = PATH + var1 + "_" + var2 + "_" + std::to_string(ctr) + ".vtk";
	std::ofstream file;
	file.open(file_name);

	file << "# vtk DataFile Version 2.0\n";
	file << "V.vtk\n";
	file << "ASCII\n";
	file << "DATASET STRUCTURED_GRID\n";
	file << "DIMENSIONS " << dim_x << " " << dim_y << " " << 1 << "\n";
	file << "Points " << dim_x * dim_y << " double\n";

	for (size_t j = 0; j < dim_y; ++j)
		for (size_t i = 0; i < dim_x; ++i)
			file << nodes_x[i] << " " << nodes_y[j] << " " << 0 << "\n";

	file << "Point_data " << vec1.size() << "\n";

	file << "\nVECTORS velocity double\n";
	for (size_t j = 0; j < dim_y; ++j) {
		for (size_t i = 0; i < dim_x; ++i)
			file << vec1[i + j * dim_x] << " " << vec2[i + j * dim_x] << " " << "0" << "\n";
		file << "\n";
	}
	file.close();
}


std::vector<double> U, V;
void interpolation() {
	U.resize((ny + 1) * (nx + 1));
	V.resize((ny + 1) * (nx + 1));
	double u_top = 0;
	// j = 0
	for (size_t i = 0; i < nx + 1; ++i) {
		if (i == 0 || i == nx)
			U[i] = 0;
		else
			U[i] = u[i] / 2.0;
		v[i] = 0;
	}

	// j from 1 to ny
	for (size_t j = 1; j < ny; ++j) {
		for (size_t i = 0; i < nx + 1; ++i) {
			if (i == 0 || i == nx)
				U[i + j * (nx + 1)] = 0;
			else
				U[i + j * (nx + 1)] = (u[i + (j - 1) * (nx + 1)] + u[i + j * (nx + 1)]) / 2.0;
		}
	}

	// j from 1 to ny
	for (size_t j = 1; j < ny; ++j) {
		for (size_t i = 0; i < nx + 1; ++i) {
			if (i == 0)
				V[i + j * (nx + 1)] = v[i + j * nx] / 2.0;
			else if (i == nx)
				V[i + j * (nx + 1)] = v[i - 1 + j * nx] / 2.0;
			else
				V[i + j * (nx + 1)] = (v[i - 1 + j * nx] + v[i + j * nx]) / 2.0;
		}
	}

	// j = ny
	for (size_t i = 0; i < nx + 1; ++i) {
		U[i + ny * (nx + 1)] = (u_top + u[i + (ny - 1) * (nx + 1)]) / 2.0;
		V[i + ny * (nx + 1)] = 0;
	}
}


int main() {
	create_mesh();


	std::string u_file_name;
	std::stringstream u_file_stream;
	size_t count_u = 0;
	std::ofstream u_file;

	std::string v_file_name;
	std::stringstream v_file_stream;
	size_t count_v = 0;
	std::ofstream v_file;

	std::string p_file_name;
	std::stringstream p_file_stream;
	size_t count_p = 0;
	std::ofstream p_file;

	std::string T_file_name;
	std::stringstream T_file_stream;
	size_t count_T = 0;
	std::ofstream T_file;

	std::string C_file_name;
	std::stringstream C_file_stream;
	size_t count_C = 0;
	std::ofstream C_file;


	int mtype = 11;												// Real unsymmetric matrix 
	// RHS and solution vectors.
	int nrhs = 1;												// Number of right hand sides. 
	// Internal solver memory pointer pt, 
	// 32-bit: int pt[64]; 64-bit: long int pt[64] 
	// or void *pt[64] should be OK on both architectures 
	void* pt[64];
	void* ptT[64];
	void* ptC[64];
	// Pardiso control parameters.
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;
	// Auxiliary variables.
	//int i;
	double ddum;												// Double dummy 
	int idum;													// Integer dummy.

	for (int i = 0; i < 64; i++) {
		iparm[i] = 0;
		pt[i] = 0;
	}

	iparm[0] = 1;												// No solver default 
	iparm[1] = 3;												// Fill-in reordering from METIS 
	// Numbers of processors, value of OMP_NUM_THREADS 
	iparm[2] = 8;
	iparm[3] = 0;												// No iterative-direct algorithm 
	iparm[4] = 0;												// No user fill-in reducing permutation 
	iparm[5] = 0;												// Write solution into x 
	iparm[6] = 0;												// Not in use 
	iparm[7] = 2;												// Max numbers of iterative refinement steps 
	iparm[8] = 0;												// Not in use 
	iparm[9] = 13;												// Perturb the pivot elements with 1E-13
	iparm[10] = 1;												// Use nonsymmetric permutation and scaling MPS 
	iparm[11] = 0;												// Not in use 
	iparm[12] = 0;												// Not in use 
	iparm[13] = 0; 												// Output: Number of perturbed pivots 
	iparm[14] = 0; 												// Not in use 
	iparm[15] = 0; 												// Not in use 
	iparm[16] = 0; 												// Not in use 
	iparm[17] = -1; 											// Output: Number of nonzeros in the factor LU 
	iparm[18] = -1; 											// Output: Mflops for LU factorization 
	iparm[19] = 0; 		                                        // Output: Numbers of CG Iterations 
	iparm[26] = 1;
	//iparm[60] = 1; //?
	maxfct = 2; 												// Maximum number of numerical factorizations. 
	mnum = 1;													// Which factorization to use. 
	//msglvl = 1;	
	msglvl = 0;												// Print statistical information in file 
	error = 0;

	/*
	iparm[0] = 1; // Нет значений по умолчанию
	iparm[1] = 2; // Используем заполнение с минимальной степенью
	iparm[7] = 1; // Максимальное количество итераций
	maxfct = 1;
	mnum = 1;
	msglvl = 1;
	error = 0;
	*/


	// Initialize error flag 
	// --------------------------------------------------------------------
	// .. Initialize the internal solver memory pointer. This is only 
	// necessary for the FIRST call of the PARDISO solver. 
	// --------------------------------------------------------------------

	for (int i = 0; i < 64; i++)
		ptT[i] = 0;

	//--------------------------------------- для T -------------------------------------------
	int mtype_T = 11;												// Real unsymmetric matrix 
	// RHS and solution vectors.
	int nrhs_T = 1;												// Number of right hand sides. 
	// Internal solver memory pointer pt, 
	// 32-bit: int pt[64]; 64-bit: long int pt[64] 
	// or void *pt[64] should be OK on both architectures 
	void* pt_T[64];
	void* ptT_T[64];
	// Pardiso control parameters.
	int iparm_T[64];
	int maxfct_T, mnum_T, phase_T, error_T, msglvl_T;
	// Auxiliary variables.
	//int i;
	double ddum_T;												// Double dummy 
	int idum_T;													// Integer dummy.

	for (int i = 0; i < 64; i++) {
		iparm_T[i] = 0;
		pt_T[i] = 0;
	}

	iparm_T[0] = 1;												// No solver default 
	iparm_T[1] = 3;												// Fill-in reordering from METIS 
	// Numbers of processors, value of OMP_NUM_THREADS 
	iparm_T[2] = 8;
	iparm_T[3] = 0;												// No iterative-direct algorithm 
	iparm_T[4] = 0;												// No user fill-in reducing permutation 
	iparm_T[5] = 0;												// Write solution into x 
	iparm_T[6] = 0;												// Not in use 
	iparm_T[7] = 2;												// Max numbers of iterative refinement steps 
	iparm_T[8] = 0;												// Not in use 
	iparm_T[9] = 13;												// Perturb the pivot elements with 1E-13
	iparm_T[10] = 1;												// Use nonsymmetric permutation and scaling MPS 
	iparm_T[11] = 0;												// Not in use 
	iparm_T[12] = 0;												// Not in use 
	iparm_T[13] = 0; 												// Output: Number of perturbed pivots 
	iparm_T[14] = 0; 												// Not in use 
	iparm_T[15] = 0; 												// Not in use 
	iparm_T[16] = 0; 												// Not in use 
	iparm_T[17] = -1; 											// Output: Number of nonzeros in the factor LU 
	iparm_T[18] = -1; 											// Output: Mflops for LU factorization 
	iparm_T[19] = 0; 		                                        // Output: Numbers of CG Iterations 
	iparm_T[26] = 1;
	//iparm[60] = 1; //?
	maxfct_T = 2; 												// Maximum number of numerical factorizations. 
	mnum_T = 1;													// Which factorization to use. 
	//msglvl = 1;	
	msglvl_T = 0;												// Print statistical information in file 
	error_T = 0;

	for (int i = 0; i < 64; i++)
		ptT_T[i] = 0;

	//--------------------------------------------------------------------------------

	for (int i = 0; i < 64; i++)
		ptC[i] = 0;

	//--------------------------------------- для C -------------------------------------------
	int mtype_C= 11;												// Real unsymmetric matrix 
	// RHS and solution vectors.
	int nrhs_C = 1;												// Number of right hand sides. 
	// Internal solver memory pointer pt, 
	// 32-bit: int pt[64]; 64-bit: long int pt[64] 
	// or void *pt[64] should be OK on both architectures 
	void* pt_C[64];
	void* ptC_C[64];
	// Pardiso control parameters.
	int iparm_C[64];
	int maxfct_C, mnum_C, phase_C, error_C, msglvl_C;
	// Auxiliary variables.
	//int i;
	double ddum_C;												// Double dummy 
	int idum_C;													// Integer dummy.

	for (int i = 0; i < 64; i++) {
		iparm_C[i] = 0;
		pt_C[i] = 0;
	}

	iparm_C[0] = 1;												// No solver default 
	iparm_C[1] = 3;												// Fill-in reordering from METIS 
	// Numbers of processors, value of OMP_NUM_THREADS 
	iparm_C[2] = 8;
	iparm_C[3] = 0;												// No iterative-direct algorithm 
	iparm_C[4] = 0;												// No user fill-in reducing permutation 
	iparm_C[5] = 0;												// Write solution into x 
	iparm_C[6] = 0;												// Not in use 
	iparm_C[7] = 2;												// Max numbers of iterative refinement steps 
	iparm_C[8] = 0;												// Not in use 
	iparm_C[9] = 13;												// Perturb the pivot elements with 1E-13
	iparm_C[10] = 1;												// Use nonsymmetric permutation and scaling MPS 
	iparm_C[11] = 0;												// Not in use 
	iparm_C[12] = 0;												// Not in use 
	iparm_C[13] = 0; 												// Output: Number of perturbed pivots 
	iparm_C[14] = 0; 												// Not in use 
	iparm_C[15] = 0; 												// Not in use 
	iparm_C[16] = 0; 												// Not in use 
	iparm_C[17] = -1; 											// Output: Number of nonzeros in the factor LU 
	iparm_C[18] = -1; 											// Output: Mflops for LU factorization 
	iparm_C[19] = 0; 		                                        // Output: Numbers of CG Iterations 
	iparm_C[26] = 1;
	//iparm[60] = 1; //?
	maxfct_C = 2; 												// Maximum number of numerical factorizations. 
	mnum_C = 1;													// Which factorization to use. 
	//msglvl = 1;	
	msglvl_C = 0;												// Print statistical information in file 
	error_C = 0;

	for (int i = 0; i < 64; i++)
		ptC_C[i] = 0;

	//--------------------------------------------------------------------------------

	size_t step = 100;

	//double mistake = 0; // разница между точным и полученным решением

	for (size_t k = 0; k < 10001; ++k) {
		matrix();
		right_part();

		// Проверка c_r
		for (size_t i = 0; i <= size_c; ++i) {
			if (c_r[i] < 1 || c_r[i] > size_m + 1) {
				std::cerr << "Ошибка в c_r[" << i << "] = " << c_r[i] << std::endl;
				return -1;
			}
		}

		// Проверка c_c
		for (size_t i = 0; i < size_m; ++i) {
			if (c_c[i] < 1 || c_c[i] > size_c) {
				std::cerr << "Ошибка в c_c[" << i << "] = " << c_c[i] << std::endl;
				return -1;
			}
		}

		int size = size_c;

		int size_T = (nx + 1) * (ny + 1);
		int size_C = size_T;

		int temp = 0;

		mnum = 1;
		phase = 11;

		if (k == 0) {
			PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, c_m, c_r, c_c, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
			if (error != 0) {
				printf("\nERROR during symbolic factorization: %d", error);
				if (error == -1) printf(" (out of memory)");
				if (error == -2) printf(" (wrong matrix indices)");
				return 1;
			}
		}
		//printf("\nReordering completed ... ");
		//printf("\nNumber of nonzeros in factors = %d", iparm[17]);
		//printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
		// --------------------------------------------------------------------
		// .. Numerical factorization.
		// --------------------------------------------------------------------
		phase = 22;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, c_m, c_r, c_c, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		if (error != 0) {
			printf("\nERROR during numerical factorization: %d", error);
			return 2;
		}
		//printf("\nFactorization completed ... ");
		// --------------------------------------------------------------------
		// .. Back substitution and iterative refinement. 
		// --------------------------------------------------------------------
		phase = 33;
		//iparm[7] = 2;
		PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &size, c_m, c_r, c_c, &idum, &nrhs, iparm, &msglvl, c_b, c_s, &error);
		if (error != 0) {
			printf("\nERROR during solution: %d", error);
			return 3;
		}

		for (size_t i = 0; i < u.size(); ++i)
			u[i] = c_s[i];
		for (size_t i = 0; i < v.size(); ++i)
			v[i] = c_s[u.size() + i];
		for (size_t i = 0; i < p.size(); ++i)
			p[i] = c_s[u.size() + v.size() + i];

		

		double div = 0.0;
		for (size_t j = 1; j < ny; ++j)
			for (size_t i = 1; i < nx; ++i) {
				//double div = (v[i + j * nx] + v[i - 1 + j * nx]) / h_x[i]
					//- (u[i + j * (nx + 1)] + u[i + (j - 1) * (nx + 1)]) / h_y[j];
				div += (u[i + 1 + j * (nx + 1)] - u[i + j * (nx + 1)]) / h_x[i]
					+ (v[i + (j + 1) * nx] - v[i + (j)*nx]) / h_y[j];
			}
		//std::cout << "\nDivergence : " << div << std::endl;


		//-----Для T-----
		matrix_T();
		right_part_T();
		mnum_T = 1;
		phase_T = 11;

		if (k == 0) {
			PARDISO(pt_T, &maxfct_T, &mnum_T, &mtype_T, &phase_T, &size_T, c_m_T, c_r_T, c_c_T,
				&idum_T, &nrhs_T, iparm_T, &msglvl_T, &ddum_T, &ddum_T, &error_T);
			if (error_T != 0) {
				printf("\nERROR during symbolic factorization: %d", error_T);
				if (error_T == -1) printf(" (out of memory)");
				if (error_T == -2) printf(" (wrong matrix indices)");
				return 1;
			}
		}

		phase_T = 22;
		PARDISO(pt_T, &maxfct_T, &mnum_T, &mtype_T, &phase_T, &size_T, c_m_T, c_r_T, c_c_T,
			&idum_T, &nrhs_T, iparm_T, &msglvl_T, &ddum_T, &ddum_T, &error_T);
		if (error_T != 0) {
			printf("\nERROR during numerical factorization: %d", error_T);
			return 2;
		}

		phase_T = 33;

		PARDISO(pt_T, &maxfct_T, &mnum_T, &mtype_T, &phase_T, &size_T, c_m_T, c_r_T, c_c_T,
			&idum_T, &nrhs_T, iparm_T, &msglvl_T, c_b_T, c_s_T, &error_T);
		if (error_T != 0) {
			printf("\nERROR during solution: %d", error_T);
			return 3;
		}

		for (size_t i = 0; i < T.size(); ++i)
			T[i] = c_s_T[i];

		//-----Для C-----
		matrix_C();
		right_part_C();
		mnum_C = 1;
		phase_C = 11;

		if (k == 0) {
			PARDISO(pt_C, &maxfct_C, &mnum_C, &mtype_C, &phase_C, &size_C, c_m_C, c_r_C, c_c_C,
				&idum_C, &nrhs_C, iparm_C, &msglvl_C, &ddum_C, &ddum_C, &error_C);
			if (error_C != 0) {
				printf("\nERROR during symbolic factorization: %d", error_C);
				if (error_C == -1) printf(" (out of memory)");
				if (error_C == -2) printf(" (wrong matrix indices)");
				return 1;
			}
		}

		phase_C = 22;
		PARDISO(pt_C, &maxfct_C, &mnum_C, &mtype_C, &phase_C, &size_C, c_m_C, c_r_C, c_c_C,
			&idum_C, &nrhs_C, iparm_C, &msglvl_C, &ddum_C, &ddum_C, &error_C);
		if (error_C != 0) {
			printf("\nERROR during numerical factorization: %d", error_C);
			return 2;
		}

		phase_C = 33;

		PARDISO(pt_C, &maxfct_C, &mnum_C, &mtype_C, &phase_C, &size_C, c_m_C, c_r_C, c_c_C,
			&idum_C, &nrhs_C, iparm_C, &msglvl_C, c_b_C, c_s_C, &error_C);
		if (error_C != 0) {
			printf("\nERROR during solution: %d", error_C);
			return 3;
		}

		for (size_t i = 0; i < C.size(); ++i)
			C[i] = c_s_C[i];



		if (k % step == 0) {
			std::cout << "\nk = " << k;
			interpolation();
			write_file("u", count_u++, u, zeros, h_y, nx + 1, ny);
			write_file("v", count_v++, v, h_x, zeros, nx, ny + 1);
			write_file("p", count_p++, p, h_x, h_y, nx, ny);
			write_vector_file("U", "V", k, U, V, nx + 1, ny + 1);

			write_file_T("T", count_T++, T, h_x, h_y, nx + 1, ny + 1);
			write_file_T("C", count_C++, C, h_x, h_y, nx + 1, ny + 1);
		}

		/*std::cout << "\n\nu:\n";
		for (size_t i = 0; i < u.size();i++)
			std::cout << u[i] << " ";

		std::cout << "\n\nv:\n";
		for (size_t i = 0; i < v.size();i++)
			std::cout << v[i] << " ";

		std::cout << "\n\np:\n";
		for (size_t i = 0; i < p.size();i++)
			std::cout << p[i] << " ";

		std::cout << "\n\nT:\n";
		for (size_t i = 0; i < T.size();i++)
			std::cout << T[i] << " ";*/


			/*
			phase = -1;
			int size_dummy = 0;
			PARDISO(
				pt, &maxfct, &mnum, &mtype, &phase,
				&size_dummy, NULL, NULL, NULL, &idum, &nrhs,
				iparm, &msglvl, NULL, NULL, &error
			);
			*/
			/*
			// Сброс указателей
			for (int i = 0; i < 64; i++) pt[i] = 0;
			*/


	}

	//std::cout << "\n Mistake is: " << mistake << "\n";

	std::cin.get();

	return 0;
}