#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <complex>
#include <valarray>


const double m1 = 2, m2 = 2;
const double l = 100.0;

const int PARAM = 4;
const int COUNT = 64;
const int STEP = 5;

double eps[6] = {0.01, -0.01, 0.02, -0.02, 0.015, -0.015};
//double eps[6] = { 0.01,0.01,0.01,0.01,0.01,0.01 };

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

class Matrix {
private:
	double* arr;
	size_t col, row;
public:
	Matrix(double** arr, int r, int c) {
		this->arr = new double[r * c];
		row = r;
		col = c;
		for (int i = 0; i < r; ++i)
			for (int j = 0; j < c; ++j)
				this->arr[i * col + j] = arr[i][j];
	}

	Matrix(size_t row, size_t col, double v = 0) {
		this->col = col;
		this->row = row;
		arr = new double[(col*row)];
		for (int i = 0; i < row * col; ++i) {
			arr[i] = 0;
		}
		for (int i = 0; i < this->getRow() && i < this->getCol(); ++i)
			arr[i * col + i] = v;
	}
	Matrix(std::vector<double> vec) {
		this->col = 1;
		this->row = vec.size();
		arr = new double[row];
		for (int i = 0; i < row; ++i)
			arr[i] = vec[i%(2*COUNT)];
	}
	void swapRows(int from, int to) {
		double tmp;
		for (int i = 0; i < col; ++i) {
			tmp = arr[from * col + i];
			arr[from * col + i] = arr[to * col + i];
			arr[to * col + i] = tmp;
		}
	}
	int getCol() const { return col; }
	int getRow() const { return row; }
	
	double operator() (int i, int j) const{
		return arr[i * col + j];
	}
	void set(size_t i, size_t j, double value) {
		arr[i * col + j] = value;
	}
	
	void add(size_t i, size_t j, double value) {
		arr[i * col + j] += value;
	}

	void print() {
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j)
				std::cout << arr[i * col + j] << ' ';
			std::cout << "\n";
		}
		std::cout << '\n';
	}
	void print(const std::string path) {
		std::ofstream file(path);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j)
				file << arr[i * col + j] << ' ';
			file << "\n";
		}
		file << '\n';
	}

	Matrix T() {
		Matrix c(this->getCol(), this->getRow());
		for (int i = 0; i < c.getRow(); ++i)
			for (int j = 0; j < c.getCol(); ++j)
				c.set(i, j, (*this)(j,i));
		return c;
	}
	Matrix& operator=(const Matrix& c) {
		if (this != &c) {
			delete[] arr;
			this->col = c.col;
			this->row = c.row;
			arr = nullptr;
			arr = new double[row * col];
			for (int i = 0; i < row * col; ++i)
				arr[i] = c.arr[i];
		}
		return *this;
	}

	Matrix operator+(const Matrix& a) const {
		Matrix c(a.getRow(), a.getCol());
		if (this->getCol() == a.getCol() && this->getRow() == a.getRow())
		{
			for (int i = 0; i < a.getRow(); ++i)
				for (int j = 0; j < a.getCol(); ++j)
					c.set(i, j, a(i,j) + (*this)(i,j));
		}
		return c;
	}

	Matrix operator-(const Matrix& a) const {
		Matrix c(a.getRow(), a.getCol());
		if (this->getCol() == a.getCol() && this->getRow() == a.getRow())
		{
			for (int i = 0; i < a.getRow(); ++i)
				for (int j = 0; j < a.getCol(); ++j)
					c.set(i, j, a(i, j) - (*this)(i, j));
		}
		return c;
	}

	Matrix operator*(const Matrix& a) const {
		Matrix c(this->getRow(), a.getCol());
		for (int i = 0; i < c.getRow(); ++i)
			for (int j = 0; j < c.getCol(); ++j)
				for (int k = 0; k < a.getRow(); ++k)
					c.add(i, j, (*this)(i, k) * a(k, j));
		return c;
	}
	Matrix operator*(const double a) const {
		Matrix c(getRow(),getCol());
		for (int i = 0; i < row * col; ++i)
			c.arr[i] = arr[i]* a;
		return c;
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
	//double dFdp[4][PARAM] =
	//{
	//{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	//{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	//{0.0,0.0,0.0,0.0,dv1dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
	//{0.0,0.0,0.0,0.0,dv2dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	//};


	//double dFdp[4][PARAM] =
	//{
	//{0.0,0.0,0.0,0.0,0.0},
	//{0.0,0.0,0.0,0.0,0.0},
	//{0.0,0.0,0.0,0.0,dv1dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
	//{0.0,0.0,0.0,0.0,dv2dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	//};
	double dFdp[4][PARAM] =
	{
	{0,0,0,0},
	{0,0,0,0},
	{0,0,0,0},
	{0,0,0,0}
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

Matrix interpolate(double t, std::vector<std::pair<double,Matrix>> deriv) {
	for (int i = 0; i < deriv.size(); ++i)
	{
		if (t > deriv[i].first)
		{
			Matrix Xt0 = deriv[i].second;
			Matrix Xt1 = deriv[i + 1].second;



			Matrix X = Xt0 + (Xt1 - Xt0) * ((t - deriv[i].first) / (deriv[i + 1].first - deriv[i].first));
			

			return X;
		}
	}
	
}


void LUP(Matrix A, Matrix& C, Matrix& P) {
	//n - размерность исходной матрицы
	const int n = A.getRow();

	C = A;

	//загружаем в матрицу P единичную матрицу
	P = Matrix(n,n,1);

	for (int i = 0; i < n; i++) {
		//поиск опорного элемента
		double pivotValue = 0;
		int pivot = -1;
		for (int row = i; row < n; row++) {
			if (fabs(C(row,i)) > pivotValue) {
				pivotValue = fabs(C(row,i));
				pivot = row;
			}
		}
		if (pivotValue != 0) {
			//меняем местами i-ю строку и строку с опорным элементом
			P.swapRows(pivot, i);
			C.swapRows(pivot, i);
			for (int j = i + 1; j < n; j++) {
				C.set(j,i, C(j,i)/C(i,i));
				for (int k = i + 1; k < n; k++)
					C.set(j,k, C(j,k)-C(j,i) * C(i,k));
			}
		}
	}
}

void Holetski(Matrix A, Matrix& L) {
	for (int i = 0; i < L.getCol(); ++i) {
		double sum = 0;
		for (int k = 0; k < i; ++k) {
			sum += std::pow(L(i, k),2);
		}
		L.set(i, i, std::sqrt(A(i, i) - sum));
		for (int j = i + 1; j < L.getRow(); ++j) {
			double s = 0;
			for (int k = 0; k < i; k++) {
				s += L(i, k) * L(j, k);
			}
			L.set(j, i, 1.0 / L(i, i) * (A(j, i) - s));
		}
	}
}

void gaussNewton(std::vector<Body> body, std::vector<std::pair<double, Matrix>> deriv, std::vector<double>& B, std::vector<double> r) {
	std::vector<double> old_B;
	for (auto a : B)
		old_B.push_back(a);
	B.clear();
	Matrix A(4 * COUNT, PARAM);
	Matrix W(4 * COUNT, 4 * COUNT,1);
	int it = 0;
	for (int i = 0; i < 4 * COUNT;)
	{
		Matrix tmp = deriv[it].second;//*/interpolate(deriv[it].first, deriv);
		for (int k = 0; k < 4; k++) {
			for (int j = 0; j < PARAM; ++j)
			{
				A.set(i, j, tmp(k,j));
			}
			i++;
		}
		it++;
	}
	//A = A * -1;
	A.print("A.txt");
	W.print("W.txt");

	Matrix At = A.T();

	Matrix AtW = At * W;
	Matrix AtWA = AtW * A;
	Matrix rB(r);
	Matrix AtWrB = AtW * rB;
	AtWA.print("AtWA.txt");
	Matrix L(PARAM, PARAM);

	//W.print();
	//
	//for (int i = 0; i < L.getCol(); ++i) {
	//	double s = AtWA(i, i);
	//	for (int k = 0; k < i; ++k)
	//		s -= std::pow(AtWA(i, k),2);
	//	L.set(i, i, std::sqrt(s));
	//	for (int j = i + 1; j < L.getRow(); ++j) {
	//		s = AtWA(j, i);
	//		for (int k = 0; k < i; ++k)
	//			s -= L(i, k) * L(j, k);
	//		L.set(j, i, s / L(i, i));
	//	}
	//}
	Holetski(AtWA, L);
	Matrix Lt = L.T();

	double* y = new double[PARAM];
	double* x = new double[PARAM];

	/// LUP 
	//Matrix P(PARAM, PARAM), C(PARAM, PARAM);
	//LUP(AtWA, C, P);
	//Matrix L(PARAM, PARAM, 1), U(PARAM,PARAM);
	//for (int i = 0; i < L.getRow(); ++i)
	//	for(int j = 0; j < i; ++j) {
	//		L.set(i, j, C(i, j));
	//}
	//for (int i = 0; i < U.getRow(); ++i)
	//	for (int j = i; j < U.getCol(); ++j) {
	//		U.set(i, j, C(i, j));
	//	}
	//Matrix LU = L * U;

	Matrix Pb = AtWrB;


	for (int i = 0; i < PARAM; ++i) {
		double sum = 0;
		for (int k = 0; k < i; ++k)
			sum += L(i, k) * y[k];
		y[i] = 1.0 / L(i, i) * (Pb(i, 0) - sum);
	}

	for (int i = PARAM - 1; i >= 0; i--) {
		double sum = 0;
		for (int k = i + 1; k < PARAM; ++k)
			sum += Lt(i, k) * x[k];
		x[i] = y[i] - sum;
	}

	for (int i = 0; i < PARAM; ++i)
		B.push_back(old_B[i] + x[i]);
	//for (int i = 0; i < 4; i++) {
	//	B[i] = old_B[i];
	//}
	//if (B[4] < 0) B[4] = old_B[4];

	//std::cout << old_B[4] << " " << B[4] << '\n';
}

void GaussNewthon(std::vector<Body> body, std::vector<std::pair<double, double**>> deriv, std::vector<double> &B, std::vector<double> r) {
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

	double** A = new double* [2 *COUNT];
	double** W = new double* [2 *COUNT];

	//double* r = new double[2 * COUNT];

	/*for (int i = 0; i < 2 * COUNT; i++) {
		r[i] = eps[rand() % 6];
	}*/

	for (int i = 0; i < 2 * COUNT; i++)
	{
		A[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			A[i][j] = 0;
	}
	for (int i = 0; i < 2 * COUNT; ++i){
		W[i] = new double[2 * COUNT];
		for (int j = 0; j < 2 * COUNT; ++j)
		{
			W[i][j] = 0;
		}
	}

	for (int i = 0; i < 2 * COUNT; ++i)
		W[i][i] = 1;
	int it = 0;
	for (int i = 0; i < 2 * COUNT;)
	{
		double** tmp = deriv[it].second;//interpolate(deriv[it].first, deriv);
		for (int k = 0; k < 2; k++) {
			for (int j = 0; j < PARAM; ++j)
			{
				A[i][j] = tmp[k][j];
			}
			i++;
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
		At[i] = new double[2 * COUNT];
		for (int j = 0; j < 2 * COUNT; ++j) {
			At[i][j] = A[j][i];
		}
	}

	double** AtW = new double* [PARAM];

	for (int i = 0; i < PARAM; i++) {
		AtW[i] = new double[2 * COUNT];
		for (int j = 0; j < 2 * COUNT; ++j)
			AtW[i][j] = 0;
	}
	for (int i = 0; i < PARAM; i++) {
		for (int j = 0; j < 2 * COUNT; j++) {
			for (int k = 0; k < 2 * COUNT; k++)
				AtW[i][j] += At[i][k] * W[k][j];
		}

	}

	double** AtWA = new double* [PARAM];
	for (int i = 0; i < PARAM; i++) {
		AtWA[i] = new double[PARAM];
		for (int j = 0; j < PARAM; j++) {
			AtWA[i][j] = 0;
			for (int k = 0; k < 2 * COUNT; ++k)
				AtWA[i][j] += AtW[i][k] * A[k][j];
		}
	}
	//double** AtA = new double* [PARAM];
	//for (int i = 0; i < PARAM; i++) {
	//	AtA[i] = new double[PARAM];
	//	for (int j = 0; j < PARAM; ++j)
	//		AtA[i][j] = 0;
	//}
	//for (int i = 0; i < PARAM; ++i)
	//	for (int j = 0; j < PARAM; j++)
	//		for (int k = 0; k < 4 * COUNT; ++k)
	//			AtA[i][j] += At[i][k] * A[k][j];
	//for (int i = 0; i < PARAM; ++i) {
	//	for (int j = 0; j < PARAM; ++j)
	//		std::cout << AtWA[i][j] << " ";
	//	std::cout << '\n';
	//}
	//std::cout << '\n';
	//double **n_AtWA = inversion(AtWA, PARAM, PARAM);

	double* rB = new double[4 * COUNT];
	for (int i = 0; i < 2 * COUNT; ++i) {
		rB[i] = r[i%(2* COUNT)];
	}

	double* AtWrB = new double[PARAM];

	for (int i = 0; i < PARAM; i++) {
		AtWrB[i] = 0;
		for (int k = 0; k < 2 * COUNT; k++) {
			AtWrB[i] += AtW[i][k] * rB[k];
		}
	}



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
		for (int j = 0; j < PARAM; ++j) {
			Lt[i][j] = L[j][i];
		}
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


	//	double** LLt = new double* [PARAM];

	//	for (int i = 0; i < PARAM; i++) {
	//		LLt[i] = new double[PARAM];
	//		for (int j = 0; j < PARAM; ++j)
	//			LLt[i][j] = 0;
	//	}
	//	for (int i = 0; i < PARAM; i++) {
	//		for (int j = 0; j < PARAM; j++) {
	//			for (int k = 0; k < PARAM; k++)
	//				LLt[i][j] += L[i][k] * Lt[k][j];
	//		}

	//	}

	//	for (int i = 0; i < PARAM; ++i) {
	//		for (int j = 0; j < PARAM; ++j)
	//			std::cout << L[i][j] << ' ';
	//		std::cout << '\t';
	//		for (int j = 0; j < PARAM; ++j)
	//			std::cout << Lt[i][j] << ' ';
	//		std::cout << '\t';
	//		for (int j = 0; j < PARAM; ++j)
	//			std::cout << LLt[i][j] << ' ';
	//		std::cout << '\n';
	//	}
	//	char a = 'a';
	//}
	std::cout << beta[4] << '\t' << B[4] << '\n';

	if (B[4] != B[4])
	{
		char a = 'a';
	}
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
	
	for (int i = 0; i < 2 * COUNT; ++i) {
		delete[] A[i], W[i];
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
	coeff _c = { k1,k2,b1,b2 };
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

	//for (int i = 0; i < 4; ++i) {
	//	delete[] km1[i], km2[i], km3[i], km4[i], stm[i];
	//}
	//delete[] km1, km2, km3, km4, stm;
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

	//for (int i = 0; i < 4; ++i) {
	//	delete[] km1[i], km2[i], km3[i], km4[i], stm[i];
	//}
	//delete[] km1, km2, km3, km4, stm;
}

bool doubleEquals(CArray arr1, CArray arr2, CArray correctArr1, CArray correctArr2) {
	double epsilon = 0.2;
	for (int i = 0; i < COUNT; i++) {
		if (std::fabs(arr1[i].real() - correctArr1[i].real()) > epsilon || std::fabs(arr2[i].real() - correctArr2[i].real()) > epsilon)
			return false;
	}

	/*for (int i = 0; i < 20; ++i)
		std::cout << correctArr1[i].real() << ' ' << arr1[i].real() << " | " << correctArr2[i].real() << ' ' << arr2[i].real() << '\n';*/
	return true;
}

int main() {
	std::vector<double> expected;
	std::vector<double> observations;
	srand(time(NULL));
	double** dXdP = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		dXdP[i] = new double[PARAM];
		for (int j = 0; j < PARAM; ++j)
			if (i == j) dXdP[i][j] = 1;
			else dXdP[i][j] = 0;
	}
	std::ofstream matr("matrix.txt");

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
	spring1.setScale(sf::Vector2f(100.0 / 502, 50.0 / 188 / 2));
	spring2.setScale(sf::Vector2f(100.0 / 502, 50.0 / 188 / 2));
	spring1.setOrigin(502 / 2, 188 / 2);
	spring2.setOrigin(502 / 2, 188 / 2);

	sf::Vertex lineX[] =
	{
		sf::Vertex(sf::Vector2f(1000,0)),
		sf::Vertex(sf::Vector2f(-1000,0))
	};
	lineX->color = sf::Color(100, 100, 100);
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

	std::vector<std::pair<double, double>> positions, correctPositions;
	bool canAdd = true;

	sf::RectangleShape obj1, obj2;
	obj1.setOrigin(25, 25);
	obj1.setSize(sf::Vector2f(50, 50));
	obj1.setFillColor(sf::Color::Red);
	obj2.setOrigin(25, 25);
	obj2.setSize(sf::Vector2f(50, 50));
	obj2.setFillColor(sf::Color::Blue);


	sf::View view = win.getDefaultView();
	view.setCenter(300, 0);
	win.setView(view);
	std::vector<std::pair<double, Matrix>> deriv;
	int step = 0;
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
		spring1.setScale((body[0].x - 25) / 502, 50.0 / 188 / 2);
		spring1.setPosition((body[0].x - 25) / 2, 0);
		spring2.setScale((body[1].x - body[0].x - 50) / 502, 50.0 / 188 / 2);
		spring2.setPosition((body[1].x + body[0].x) / 2, 0);

		win.draw(obj1);
		win.draw(spring1);
		win.draw(obj2);
		win.draw(spring2);

		RungeKutta(body, dXdP, h, f);
		step++;
		t += h;
		if (step % STEP == 0)
		{
			expected.push_back(body[0].x);
			expected.push_back(body[1].x);
			observations.push_back(body[0].x);
			observations.push_back(body[1].x);
			file << t << " " << body[0].x << " " << body[1].x << '\n';
			positions.push_back(std::make_pair(body[0].x, body[1].x));
			if (canAdd)	correctPositions.push_back(std::make_pair(body[0].x, body[1].x));
			Matrix tmp(dXdP, 4, PARAM);
			deriv.push_back(std::make_pair(t, tmp));
			//for (int i = 0; i < 4; ++i) {
			//	for (int j = 0; j < PARAM; ++j)
			//		matr << tmp[i][j] << " ";
			//	matr << '\n';
			//}
			//matr << '\n';
		}

		win.display();

		if (step / STEP >= COUNT)
			win.close();
	}

	file.close();
	canAdd = false;
	std::vector<Body> bodyT;
	Complex correctArr1[COUNT], correctArr2[COUNT];
	for (int i = 0; i < COUNT; i++)
	{
		correctArr1[i] = correctPositions[i].first;
		correctArr2[i] = correctPositions[i].second;
	}
	Complex arr1[COUNT], arr2[COUNT];
	CArray data1, data2, correctData1, correctData2;
	std::vector<double> B = { 145,355,2,3,6,6,0.05,0.08 };
	int k = 0;
	bodyT = body;
	std::vector<std::pair<double, Matrix>> derivTmp;
	derivTmp = deriv;
	std::vector<double> r;
	while (true) {
		//for (int p = 0; p < derivTmp.size(); ++p) {
		//	derivTmp[p].second = derivTmp[p].second*0.01;
		//}

		for (int i = 0; i < expected.size(); ++i)
			r.push_back(expected[i] - observations[i] - eps[rand() % 6]);
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

		//for (int k = 0; k < 20; ++k) {
		//	for (int i = 0; i < 4; ++i) {
		//		for (int j = 0; j < PARAM; ++j)
		//			std::cout << derivTmp[k].second(i,j) << ' ';
		//		std::cout << '\n';
		//	}
		//	std::cout << '\n';
		//}

		gaussNewton(bodyT, derivTmp, B, r);

		//for (auto a : B)
		//	std::cout << a << ' ';
		//std::cout << "\n\n";

		if (k % 1 == 0) {
			std::cout << "Debug: ";
			for (int i = 0; i < PARAM; ++i)
				std::cout << B[i] << ' ';
			std::cout << '\n';
		}
		for (int i = 0; i < COUNT; i++)
		{
			arr1[i] = positions[i].first;
			arr2[i] = positions[i].second;
		}
		data1 = CArray(arr1, COUNT);
		data2 = CArray(arr2, COUNT);
		correctData1 = CArray(correctArr1, COUNT);
		correctData2 = CArray(correctArr2, COUNT);
		//fft(data1);
		//fft(data2);
		//fft(correctData1);
		//fft(correctData2);
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
		//	std::cout << "Debug: ";
		//	for (int j = 0; j < 5; j++)
		//		std::cout << dXdPT[i][j] << " ";
		//	std::cout << "\n";
		//}
		//std::cout << "----------\n";
		bodyT.clear();
		bodyT.push_back(Body(B[0], B[2]));
		bodyT.push_back(Body(B[1], B[3]));
		derivTmp.clear();
		positions.clear();
		r.clear();
		observations.clear();

		//std::cout << k << "###############################\n";
		int st = 0;
		while (st / STEP <= COUNT) {
			st++;
			RungeKuttaT(bodyT, dXdPT, h, fT, B);


			//for (int i = 0; i < 4; ++i) {
			//	for (int j = 0; j < 5; j++)
			//		std::cout << dXdPT[i][j] << " ";
			//	std::cout << "\n";
			//}
			//std::cout << '\n';
			ta += h;

			if (st % STEP == 0) {
				observations.push_back(bodyT[0].x);
				observations.push_back(bodyT[1].x);
				positions.push_back(std::make_pair(bodyT[0].x, bodyT[1].x));
				Matrix tmp(dXdPT, 4, PARAM);
				derivTmp.push_back(std::make_pair(ta, tmp));
			}

		}



		char a = 'a';
		//if (dXdPT[0][0] != dXdPT[0][0]) std::cout << "First\n";
	}
	for (int i = 0; i < PARAM; ++i)
		std::cout << B[i] << ' ';
	std::cout << '\n';
	//std::cout << k << '\n';
	return 0;
}
