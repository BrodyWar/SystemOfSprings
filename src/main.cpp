#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <complex>
#include <valarray>


const double m1 = 2, m2 = 2;
const double l = 100.0;

const int PARAM = 5;
const int COUNT = 20;

double eps[6] = {0.01, -0.01, 0.02, -0.02, 0.015, -0.015};

struct coeff {
	double k1, k2, b1, b2;
};

class Body {
public:
	double x, v, a;
	Body() {
		x = 0;
		v = 0;
		a = 0;
	}
	Body(double x, double v) :x(x), v(v), a(0){}
	Body operator+(const Body& a)const {
		return Body{ this->x + a.x, this->v + a.v };
	}
	Body operator-(const Body& a)const {
		return Body{ this->x - a.x, this->v - a.v };
	}
	Body operator*(const Body& a)const {
		return Body{ this->x * a.x, this->v * a.v };
	}
	Body operator*(const double a)const {
		return Body{ this->x * a, this->v * a };
	}
	Body operator/(const double a)const {
		return Body{ this->x / a, this->v / a };
	}
};

double dv1dx1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return (-k1 / m1 - k2 / m1);
}
double dv1dx2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return k2 / m1;
}
double dv1dv1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -b1 / m1;
}
double dv1dv2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return b2 / m1;
}
double dv1dk1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return (l - x1) / m1;
}
double dv1dk2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -(l + x1 - x2) / m1;
}
double dv1db1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -v1 / m1;
}
double dv1db2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return v2 / m1;
}
double dv2dx1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return k2 / m2;
}
double dv2dx2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -k2 / m2;
}
double dv2dv1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return 0;
}
double dv2dv2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -b2 / m2;
}
double dv2dk1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return 0;
}
double dv2dk2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return (l + x1 - x2) / m2;
}
double dv2db1(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return 0;
}
double dv2db2(double x1, double x2, double v1, double v2, double k1, double k2, double b1, double b2) {
	return -v2 / m2;
}


const double PI = 3.141592653589793238460;

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

void fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N / 2, 2)];
	CArray  odd = x[std::slice(1, N / 2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t k = 0; k < N / 2; ++k)
	{
		Complex t = std::polar(1.0, -2 * PI * k / N) * odd[k];
		x[k] = even[k] + t;
		x[k + N / 2] = even[k] - t;
	}
}


void updateSystem(std::vector<Body> body, double** dxdp, double** ndxdp, coeff _c)
{
	/*double dFdp[4][8] =
	{
	{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,dv1dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
	{0.0,0.0,0.0,0.0,dv2dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	};*/


	double dFdp[4][PARAM] =
	{
	{0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,dv1dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
	{0.0,0.0,0.0,0.0,dv2dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	};

	double dFdX[4][4] =
	{
		{0,0,1,0},
		{0,0,0,1},
		{dv1dx1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dx2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dv1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dv2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
		{dv2dx1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dx2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dv1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dv2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	};
	double** c = new double* [4];
	for (int i = 0; i < 4; i++) {
		c[i] = new double[PARAM];
		for (int j = 0; j < PARAM; j++) {
			c[i][j] = 0;
			for (int k = 0; k < 4; k++)
				c[i][j] += dFdX[i][k] * dxdp[k][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < PARAM; j++) {
			ndxdp[i][j] = dFdp[i][j] + c[i][j];
		}
	}

	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 5; ++j)
	//		std::cout << dFdp[i][j] << ' ';
	//	std::cout << "\t\t";
	//	for (int j = 0; j < 5; ++j)
	//		std::cout << dxdp[i][j] << ' ';
	//	std::cout << '\n';
	//}
	//std::cout << _c.k1 << '\n';
	//std::cout << '\n';

}

double** interpolate(double t, std::vector<std::pair<double,double**>> deriv) {
	for (int i = 0; i < deriv.size()-1; ++i)
	{
		if (std::fabs(deriv[i].first - t) <= 1.2)
		{
			double** Xt0 = deriv[i].second;
			double** Xt1 = deriv[i + 1].second;

			double** X = new double* [4];
			for (int k = 0; k < 4; k++) {
				X[k] = new double[PARAM];
				for (int h = 0; h < PARAM; h++)
					X[k][h] = 0;
			}


			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < PARAM; k++)
					X[j][k] = Xt0[j][k] + (Xt1[j][k] - Xt0[j][k]) * ((t - deriv[i].first) / (deriv[i + 1].first - deriv[i].first));
			}

			return X;
		}
	}
}



void GaussNewthon(std::vector<Body> body, std::vector<std::pair<double, double**>> deriv, std::vector<double> &B) {
	double* beta = new double[PARAM];

	//if (deriv[0].second[0][0] != deriv[0].second[0][0]) std::cout << "Second!!!\n";
	{
		for (int i = 0; i < PARAM; i++)
		beta[2] = B[2];
		beta[3] = B[3];
		beta[0] = B[0];
		beta[1] = B[1];
		beta[4] = B[4];
		//beta[5] = B[5];
		//beta[6] = B[6];
		//beta[7] = B[7];
	}
	B.clear();

	double** A = new double* [COUNT];
	double** W = new double* [COUNT];

	double* r = new double[COUNT];

	for (int i = 0; i < COUNT; i++) {
		r[i] = eps[rand() % 6];
	}

	for (int i = 0; i < COUNT; i++)
	{
		A[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			A[i][j] = 0;
	}
	for (int i = 0; i < COUNT; ++i){
		W[i] = new double[COUNT];
		for (int j = 0; j <  COUNT; ++j)
		{
			W[i][j] = 0;
		}
	}

	for (int i = 0; i < COUNT; ++i)
		W[i][i] = 1;
	int it = 0;
	for (int i = 0; i <COUNT; i++)
	{
		double** tmp = interpolate(deriv[it].first, deriv);
		for (int j = 0; j < PARAM; ++j)
		{
			A[i][j] = tmp[i%4][j];
		}
		it++;

		//if (tmp[0][0] != tmp[0][0]) {
		//	for (int k = 0; k < 20; ++k) {
		//		for (int i = 0; i < 4; ++i) {
		//			for (int j = 0; j < 5; ++j)
		//				std::cout << deriv[k].second[i][j] << " ";
		//			std::cout << '\n';
		//		}
		//		std::cout << '\n';
		//	}
		//}
	}
	double** At = new double* [PARAM];
	for (int i = 0; i < PARAM; i++)
	{
		At[i] = new double[COUNT];
		for (int j = 0; j < COUNT; ++j)
			At[i][j] = A[j][i];
	}

	double** AtW = new double* [PARAM];

	for (int i = 0; i < PARAM; i++) {
		AtW[i] = new double[ COUNT];
		for (int j = 0; j <  COUNT; ++j)
			AtW[i][j] = 0;
	}
	for (int i = 0; i < PARAM; i++) {
		for (int j = 0; j < COUNT; j++) {
			for (int k = 0; k < COUNT; k++)
				AtW[i][j] += At[i][k] * W[k][j];
		}

	}

	double** AtWA = new double* [PARAM];
	for (int i = 0; i < PARAM; i++) {
		AtWA[i] = new double[PARAM];
		for (int j = 0; j < PARAM; j++) {
			AtWA[i][j] = 0;
			for (int k = 0; k < COUNT; ++k)
				AtWA[i][j] += AtW[i][k] * A[k][j];
		}
	}
	//double **n_AtWA = inversion(AtWA, PARAM, PARAM);

	double** rB = new double* [COUNT];
	for (int i = 0; i < COUNT; ++i) {
		rB[i] = new double[1];
		rB[i][0] = r[i];
	}

	double* AtWrB = new double[PARAM];

	for (int i = 0; i < PARAM; i++) {
		AtWrB[i] = 0;
		for (int k = 0; k < COUNT; k++) {
			AtWrB[i] += AtW[i][k] * rB[k][0];
		}
	}
	//for (int i = 0; i < PARAM; ++i) {
	//	for (int j = 0; j < PARAM; ++j)
	//		std::cout << AtWA[i][j] << ' ';
	//	std::cout << '\n';
	//}
	//std::cout << '\n';


	double** L = new double* [PARAM];
	double** Lt = new double* [PARAM];
	for (int i = 0; i < PARAM; ++i) {
		L[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			L[i][j] = 0;
	}
	L[0][0] = std::sqrt(AtWA[0][0]);
	for (int j = 1; j < PARAM; ++j)
		L[j][0] = AtWA[j][0] / L[0][0];
	
	L[1][1] = std::sqrt(AtWA[1][1] - std::pow(L[1][0], 2));
	for (int i = 1; i < PARAM - 1; ++i) {
		{
			double sum = 0;
			for (int p = 0; p <= i - 1; ++p) {
				sum += std::pow(L[i][p], 2);
			}
			L[i][i] = std::sqrt(AtWA[i][i] - sum);
		}

		for (int j = i + 1; j < PARAM; ++j) {
			double sum = 0;
			for (int p = 0; p <= i - 1; ++p) {
				sum += L[i][p] * L[j][p];
			}
			L[j][i] = 1.0 / L[i][i] * (AtWA[j][i] - sum);
		}
	}
	double sum = 0;
	for (int p = 0; p < PARAM-1; ++p)
		sum += std::pow(L[PARAM-1][p], 2);
	L[PARAM - 1][PARAM - 1] = std::sqrt(AtWA[PARAM - 1][PARAM - 1] - sum);


	for (int i = 0; i < PARAM; i++)
	{
		Lt[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			Lt[i][j] = L[j][i];
	}

	double* y = new double[PARAM];
	double* x = new double[PARAM];

	for (int i = 0; i < PARAM; ++i) {
		double sum = 0;
		for (int k = 0; k < i; ++k)
			sum += L[i][k] * y[k];
		y[i] = 1.0 / L[i][i] * (AtWrB[i] - sum);
	}

	for (int i = PARAM-1; i >=0; i--) {
		double sum = 0;
		for (int k = i + 1; k < PARAM; ++k)
			sum += Lt[i][k] * x[k];
		x[i] = y[i] - sum;
	}

	for (int i = 0; i < PARAM; ++i)
		B.push_back(beta[i] + x[i]);


	if (B[0] != B[0])
	{
		char a = 'a';
	}

	//double** LLt = new double* [PARAM];

	//for (int i = 0; i < PARAM; i++) {
	//	LLt[i] = new double[PARAM];
	//	for (int j = 0; j < PARAM; ++j)
	//		LLt[i][j] = 0;
	//}
	//for (int i = 0; i < PARAM; i++) {
	//	for (int j = 0; j < PARAM; j++) {
	//		for (int k = 0; k < PARAM; k++)
	//			LLt[i][j] += L[i][k] * Lt[k][j];
	//	}

	//}

	//for (int i = 0; i < PARAM; ++i) {
	//	for (int j = 0; j < PARAM; ++j)
	//		std::cout << L[i][j] << ' ';
	//	std::cout << '\t';
	//	for (int j = 0; j < PARAM; ++j)
	//		std::cout << Lt[i][j] << ' ';
	//	std::cout << '\t';
	//	for (int j = 0; j < PARAM; ++j)
	//		std::cout << LLt[i][j] << ' ';
	//	std::cout << '\n';
	//}
	//std::cout << '\n';
	/*double** mTmp = new double* [PARAM];
	for (int i = 0; i < PARAM; i++) {
		mTmp[i] = new double[1];
		mTmp[i][0] = 0;
		for (int k = 0; k < PARAM; k++)
			mTmp[i][0] += n_AtWA[i][k] * AtWrB[k][0];
	}*/

	//double** rBeta = new double* [PARAM];
	//for (int i = 0; i < PARAM; i++) {
	//	rBeta[i] = new double[1];
	//	rBeta[i][0] = beta[i][0] - mTmp[i][0];
	//	B.push_back(rBeta[i][0]);
	//}

	std::cout << beta[0] << '\t' << B[0] << '\n';

	//if (B[0] != B[0]) {
	//	for (int i = 0; i < 40; ++i) {
	//		for (int j = 0; j < 1; ++j)
	//			std::cout << rB[i][j] << " ";
	//		std::cout << '\n';
	//	}
	//	std::cout << '\n';	
	//}


	//if (B[0] != B[0]) {
	//	for (int k = 0; k < 20; ++k) {
	//		for (int i = 0; i < 4; ++i) {
	//			for (int j = 0; j < 5; ++j)
	//				std::cout << deriv[k].second[i][j] << ' ';
	//			std::cout << '\n';
	//		}
	//		std::cout << '\n';
	//	}
	//}
	
	for (int i = 0; i < COUNT; ++i) {
		delete[] A[i], W[i], r[i], rB[i];
	}
	delete[] A, W, r, rB;

	for (int i = 0; i < PARAM; ++i) {
		delete[] At[i], AtW[i], AtWA[i], AtWrB[i];
	}
	delete[] At, AtW, AtWA, AtWrB;
}

void f(std::vector<Body>& body, std::vector<Body>& dot, double** ndxdp, double** dxdp) {
	double k1 = 6, k2 = 6, b1 = 0.05, b2 = 0.08;

	dot[0].x = body[0].v;
	dot[0].v = (-k1 * (body[0].x - l) - b1 * body[0].v + k2 * (body[1].x - body[0].x - l) + b2 * body[1].v) / m1;
	dot[1].x = body[1].v;
	dot[1].v = (-k2 * (body[1].x - body[0].x - l) - b2 * body[1].v) / m2;
	coeff _c = {k1,k2,b1,b2};
	updateSystem(body, dxdp, ndxdp, _c);
}

void fT(std::vector<Body>& body, std::vector<Body>& dot, double** ndxdp, double** dxdp, std::vector<double> B) {
	double k1 = 6, k2 = 6, b1 = 0.05, b2 = 0.08;
	//coeff _c = { B[4],B[5],B[6],B[7] };
	coeff _c = { B[4],k2,b1,b2 };
	dot[0].x = body[0].v;
	dot[0].v = (-_c.k1 * (body[0].x - l) - _c.b1 * body[0].v + _c.k2 * (body[1].x - body[0].x - l) + _c.b2 * body[1].v) / m1;
	dot[1].x = body[1].v;
	dot[1].v = (-_c.k2 * (body[1].x - body[0].x - l) - _c.b2 * body[1].v) / m2;
	
	updateSystem(body, dxdp, ndxdp, _c);
}

void step(std::vector<Body>& st, std::vector<Body>& body, std::vector<Body> k,   double** dxdp, double** km, double** ndxdp, double h) {
	for (int i = 0; i < st.size(); i++) {
		st[i] = body[i] + k[i] * h / 2;
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < PARAM; ++j)
			ndxdp[i][j] = dxdp[i][j] + h / 2 * km[i][j];
}



void RungeKutta(std::vector<Body>& body, double** dxdp, double h, void (*f)(std::vector<Body>&, std::vector<Body>&, double**, double**)) {
	std::vector<Body> k1, k2, k3, k4, st;
	double **km1=new double*[4],
		   **km2=new double*[4], 
		   **km3=new double*[4], 
		   **km4=new double*[4], 
		   **stm=new double*[4];
	for (int i = 0; i < 4; ++i)
	{
		km1[i] = new double[PARAM],
		km2[i] = new double[PARAM],
		km3[i] = new double[PARAM],
		km4[i] = new double[PARAM],
		stm[i] = new double[PARAM];

		for (int j = 0; j < PARAM; ++j)
		{
			km1[i][j] = 0;
			km2[i][j] = 0;
			km3[i][j] = 0;
			km4[i][j] = 0;
			stm[i][j] = 0;
		}
	}

	for (int i = 0; i < body.size(); ++i)
	{
		k1.push_back(Body());
		k2.push_back(Body());
		k3.push_back(Body());
		k4.push_back(Body());
		st.push_back(Body());
	}
	f(body,k1,km1,dxdp);
	step(st,body,k1,dxdp,km1,stm,h);
	f(st,k2,km2,stm);

	step(st,body,k2,dxdp,km2,stm,h);
	f(st,k3,km3,stm);

	step(st,body,k3,dxdp,km3,stm,h*2);
	f(st,k4,km4,stm);

	for (int i = 0; i < body.size(); ++i) {
		body[i] = body[i] + ((k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) * h / 6);
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < PARAM; ++j)
			dxdp[i][j] += ((km1[i][j] + km2[i][j] * 2 + km3[i][j] * 2 + km4[i][j]) * h / 6);

	for (int i = 0; i < 4; ++i) {
		delete[] km1[i], km2[i], km3[i], km4[i], stm[i];
	}
	delete[] km1, km2, km3, km4, stm;
}

void RungeKuttaT(std::vector<Body>& body, double** dxdp, double h, void (*fT)(std::vector<Body>&, std::vector<Body>&, double**, double**, std::vector<double> B), std::vector<double> B) {
	std::vector<Body> k1, k2, k3, k4, st;
	double** km1 = new double* [4],
		** km2 = new double* [4],
		** km3 = new double* [4],
		** km4 = new double* [4],
		** stm = new double* [4];
	for (int i = 0; i < 4; ++i)
	{
			km1[i] = new double[5],
			km2[i] = new double[5],
			km3[i] = new double[5],
			km4[i] = new double[5],
			stm[i] = new double[5];

		for (int j = 0; j < PARAM; ++j)
		{
			km1[i][j] = 0;
			km2[i][j] = 0;
			km3[i][j] = 0;
			km4[i][j] = 0;
			stm[i][j] = 0;
		}
	}

	for (int i = 0; i < body.size(); ++i)
	{
		k1.push_back(Body());
		k2.push_back(Body());
		k3.push_back(Body());
		k4.push_back(Body());
		st.push_back(Body());
	}
	fT(body, k1, km1, dxdp,B);
	step(st, body, k1, dxdp, km1, stm, h);
	fT(st, k2, km2, stm, B);
	//std::cout << "1----------\n";

	step(st, body, k2, dxdp, km2, stm, h);
	fT(st, k3, km3, stm, B);
	//std::cout << "2----------\n";
	step(st, body, k3, dxdp, km3, stm, h * 2);
	fT(st, k4, km4, stm, B);
	//std::cout << "3----------\n";
	for (int i = 0; i < body.size(); ++i) {
		body[i] = body[i] + ((k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) * h / 6);
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < PARAM; ++j)
			dxdp[i][j] += ((km1[i][j] + km2[i][j] * 2 + km3[i][j] * 2 + km4[i][j]) * h / 6);

	for (int i = 0; i < 4; ++i) {
		delete[] km1[i], km2[i], km3[i], km4[i], stm[i];
	}
	delete[] km1, km2, km3, km4, stm;
}

bool doubleEquals(CArray arr1, CArray arr2, CArray correctArr1, CArray correctArr2) {
	double epsilon = 100;
	for (int i = 0; i < 20; i++) {
		if (std::fabs(arr1[i].real() - correctArr1[i].real()) > epsilon || std::fabs(arr2[i].real() - correctArr2[i].real()) > epsilon)
			return false;
	}

	/*for (int i = 0; i < 20; ++i)
		std::cout << correctArr1[i].real() << ' ' << arr1[i].real() << " | " << correctArr2[i].real() << ' ' << arr2[i].real() << '\n';*/
	return true;
}

int main() {
	srand(time(NULL));
	double** dXdP = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		dXdP[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			if (i == j) dXdP[i][j] = 1;
			else dXdP[i][j] = 0;
	}
	std::ofstream matrix("matrix.txt");

	sf::RenderWindow win(sf::VideoMode(980, 720), "Window");
	win.setFramerateLimit(60);
	std::ofstream file("time.txt");

	sf::Image spring_image;
	spring_image.loadFromFile("images/spring.png");
	sf::Texture spring_texture;
	spring_texture.loadFromImage(spring_image);
	sf::Sprite spring1, spring2;
	spring1.setTexture(spring_texture);
	spring2.setTexture(spring_texture);
	spring1.setTextureRect(sf::IntRect(0, 0, 502, 188));
	spring2.setTextureRect(sf::IntRect(0, 0, 502, 188));
	spring1.setScale(sf::Vector2f(100.0/502,50.0/188/2));
	spring2.setScale(sf::Vector2f(100.0 / 502, 50.0 / 188/2));
	spring1.setOrigin(502/2, 188/2);
	spring2.setOrigin(502/2, 188/2);

	sf::Vertex lineX[] =
	{
		sf::Vertex(sf::Vector2f(1000,0)),
		sf::Vertex(sf::Vector2f(-1000,0))
	};
	lineX->color= sf::Color(100, 100, 100);
	sf::Vertex lineY[] =
	{
		sf::Vertex(sf::Vector2f(0,-1000)),
		sf::Vertex(sf::Vector2f(0,1000))
	};
	lineY->color = sf::Color(100, 100, 100);

	double h = 0.1;
	double t = 0;

	std::vector<Body> body;
	body.push_back(Body(150, 0));
	body.push_back(Body(350, 0));

	std::vector<std::pair<double,double>> positions, correctPositions;
	bool canAdd = true;

	sf::RectangleShape obj1, obj2;
	obj1.setOrigin(25,25);
	obj1.setSize(sf::Vector2f(50, 50));
	obj1.setFillColor(sf::Color::Red);
	obj2.setOrigin(25, 25);
	obj2.setSize(sf::Vector2f(50, 50));
	obj2.setFillColor(sf::Color::Blue);


	sf::View view = win.getDefaultView();
	view.setCenter(300, 0);
	win.setView(view);
	std::vector<std::pair<double, double**>> deriv;
	while (win.isOpen()) {
		sf::Event e;
		while (win.pollEvent(e)) {
			if (e.type == sf::Event::Closed)
				win.close();
		}
		win.clear();

		win.draw(lineX, 2, sf::Lines);
		win.draw(lineY, 2, sf::Lines);

		obj1.setPosition(body[0].x, 0);
		obj2.setPosition(body[1].x, 0);
		spring1.setScale((body[0].x - 25) / 502, 50.0 / 188/2);
		spring1.setPosition((body[0].x - 25)/2, 0);
		spring2.setScale((body[1].x - body[0].x -50) / 502, 50.0 / 188/ 2);
		spring2.setPosition((body[1].x + body[0].x) / 2, 0);

		win.draw(obj1);
		win.draw(spring1);
		win.draw(obj2);
		win.draw(spring2);

		RungeKutta(body, dXdP, h, f);
		t += h;
		if (t - int(t) < h)
		{
			file << t << " " << body[0].x << " " << body[1].x << '\n';
			positions.push_back(std::make_pair(body[0].x, body[1].x));
			if (canAdd)	correctPositions.push_back(std::make_pair(body[0].x, body[1].x));
			deriv.push_back(std::make_pair(t, dXdP));
			for (int i = 0; i < 4; ++i) {
				for (int j = 0; j < PARAM; ++j)
					matrix << dXdP[i][j] << " ";
				matrix << '\n';
			}
			matrix << '\n';
		}
		
		win.display();

		if (t >= 20)
			win.close();
	}

	file.close();
	canAdd = false;
	std::vector<Body> bodyT;
	Complex correctArr1[20], correctArr2[20];
	for (int i = 0; i < 20; i++)
	{
		correctArr1[i] = correctPositions[i].first;
		correctArr2[i] = correctPositions[i].second;
	}
	Complex arr1[20], arr2[20];
	CArray data1, data2, correctData1, correctData2;
	std::vector<double> B = {150,350,0,0,5.8,6,0.05,0.08};
	int k = 0;
	bodyT = body;
	std::vector<std::pair<double, double**>> derivTmp;
	derivTmp = deriv;
	while (true) {

		//for (int k = 0; k < 20; k++) {
		//	for (int i = 0; i < 4; ++i) {
		//		for (int j = 0; j < 5; j++)
		//			std::cout << derivTmp[k].second[i][j] << " ";
		//		std::cout << "\n";
		//	}
		//	std::cout << '\n';
		//}
		//std::cout << '\n';

		//for (auto a : B)
		//	std::cout << a << ' ';
		//std::cout << "\n---------------------\n";

		GaussNewthon(bodyT, derivTmp, B);

		//for (auto a : B)
		//	std::cout << a << ' ';
		//std::cout << "\n\n";

		 //if (k % 1 == 0) std::cout << "\\" << B[4] /*<< ' ' << B[5] << ' ' << B[6] << ' ' << B[7]*/ << '\n';
		for (int i = 0; i < 20; i++)
		{
			arr1[i] = positions[i].first;
			arr2[i] = positions[i].second;
		}
		data1 = CArray(arr1, 20);
		data2 = CArray(arr2, 20);
		correctData1 = CArray(correctArr1, 20);
		correctData2 = CArray(correctArr2, 20);
		fft(data1);
		fft(data2);
		fft(correctData1);
		fft(correctData2);
		//std::cout << data1[0].real() << " " << data2[0].real() << " " << correctData1[0].real() << " " << correctData2[0].real() << '\n';
		if (k != 0 && doubleEquals(data1, data2, correctData1, correctData2))
			break;
		double ta = 0;
		k++;
		double** dXdPT = new double* [4];
		for (int i = 0; i < 4; i++)
		{
			dXdPT[i] = new double[PARAM];
			for (int j = 0; j < PARAM; ++j)
				if (i == j) dXdPT[i][j] = 1;
				else dXdPT[i][j] = 0;
		}

		//for (int i = 0; i < 4; ++i) {
		//	for (int j = 0; j < 5; j++)
		//		std::cout << dXdPT[i][j] << " ";
		//	std::cout << "\n";
		//}
		//std::cout << "----------\n";
		bodyT.clear();
		bodyT.push_back(Body(150, 0));
		bodyT.push_back(Body(350, 0));
		derivTmp.clear();
		positions.clear();

		//std::cout << k << "###############################\n";
		while (ta <= 20) {
			RungeKuttaT(bodyT, dXdPT, h, fT,B);


			//for (int i = 0; i < 4; ++i) {
			//	for (int j = 0; j < 5; j++)
			//		std::cout << dXdPT[i][j] << " ";
			//	std::cout << "\n";
			//}
			//std::cout << '\n';
			ta += h;

			if (ta - (int)ta > h) {
				positions.push_back(std::make_pair(bodyT[0].x, bodyT[1].x));
				derivTmp.push_back(std::make_pair(ta, dXdPT));
			}
		}
		//if (dXdPT[0][0] != dXdPT[0][0]) std::cout << "First\n";
	}
	std::cout << B[4] /*<< ' ' << B[5] << ' ' << B[6] << ' ' << B[7]*/ << '\n';
	//std::cout << k << '\n';
	return 0;
}
