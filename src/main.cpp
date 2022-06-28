#include <iostream>
#include <SFML/Graphics.hpp>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <complex>
#include <valarray>


double m1 = 2, m2 = 2;
double l = 100.0;

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
	double dFdp[4][8] =
	{
	{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
	{0.0,0.0,0.0,0.0,dv1dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
	{0.0,0.0,0.0,0.0,dv2dk1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dk2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2db2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	};

	double dFdX[4][4] =
	{
		{0.0,0.0,1.0,0.0},
		{0.0,0.0,0.0,1.0},
		{dv1dx1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dx2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dv1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv1dv2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)},
		{dv2dx1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dx2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dv1(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2),dv2dv2(body[0].x, body[1].x, body[0].v, body[1].v, _c.k1, _c.k2, _c.b1, _c.b2)}
	};
	double** c = new double* [4];
	for (int i = 0; i < 4; i++) {
		c[i] = new double[8];
		for (int j = 0; j < 8; j++) {
			c[i][j] = 0;
			for (int k = 0; k < 4; k++)
				c[i][j] += dFdX[i][k] * dxdp[k][j];
		}
	}
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 8; j++) {
			ndxdp[i][j] = dFdp[i][j] + c[i][j];
		}
	}

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
				X[k] = new double[8];
				for (int h = 0; h < 8; h++)
					X[k][h] = 0;
			}


			for (int j = 0; j < 4; j++)
			{
				for (int k = 0; k < 8; k++)
					X[j][k] = Xt0[j][k] + (Xt1[j][k] - Xt0[j][k]) * ((t - deriv[i].first) / (deriv[i + 1].first - deriv[i].first));
			}
			return X;
		}
	}
}

double** inversion(double** A, int n, int m)
{
	double temp;

	double** E = new double*[n];
	for (int i = 0; i < n; i++)
	{
		E[i] = new double[m];
		for (int j = 0; j < m; ++j)
			E[i][j] = 0;
	}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < n; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < m; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < m; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < m; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = m - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < n; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			A[i][j] = E[i][j];

	return A;
}


void GaussNewthon(std::vector<Body> body, std::vector<std::pair<double, double**>> deriv, std::vector<double> &B) {
	double** beta = new double* [8];
	{
		for (int i = 0; i < 8; i++)
			beta[i] = new double[1];
		beta[2][0] = B[2];
		beta[3][0] = B[3];
		beta[0][0] = B[0];
		beta[1][0] = B[1];
		beta[4][0] = B[4];
		beta[5][0] = B[5];
		beta[6][0] = B[6];
		beta[7][0] = B[7];
	}
	B.clear();

	double** A = new double* [40];
	double** W = new double* [40];

	double* r = new double[40];

	for (int i = 0; i < 40; i++) {
		r[i] = eps[rand() % 6];
	}

	for (int i = 0; i < 40; i++)
	{
		A[i] = new double[8];
		for (int j = 0; j < 8; ++j)
			A[i][j] = 0;
	}
	for (int i = 0; i < 40; ++i){
		W[i] = new double[40];
		for (int j = 0; j < 40; ++j)
		{
			W[i][j] = 0;
			if (i == j) W[i][j] = 1;
		}
	}
	int it = 0;
	for (int i = 0; i < 40; i++)
	{
		double** tmp = interpolate(deriv[it%20].first, deriv);
		for (int j = 0; j < 8; ++j)
		{
			A[i][j] = tmp[i%4][j];
		}
		it++;
	}
	double** At = new double* [8];
	for (int i = 0; i < 8; i++)
	{
		At[i] = new double[40];
		for (int j = 0; j < 40; ++j)
			At[i][j] = A[j][i];
	}

	double** AtW = new double* [8];
	for (int i = 0; i < 8; i++) {
		AtW[i] = new double[40];
		for (int j = 0; j < 40; j++) {
			AtW[i][j] = 0;
			for (int k = 0; k < 40; k++)
				AtW[i][j] += At[i][k] * W[k][j];
		}
	}

	double** AtWA = new double* [8];
	for (int i = 0; i < 8; i++) {
		AtWA[i] = new double[8];
		for (int j = 0; j < 8; j++) {
			AtWA[i][j] = 0;
			for (int k = 0; k < 40; ++k)
				AtWA[i][j] += AtW[i][k] * A[k][j];
		}
	}

	AtWA = inversion(AtWA, 8, 8);

	double** rB = new double* [40];
	for (int i = 0; i < 40; ++i) {
		rB[i] = new double[1];
		rB[i][0] = r[i];
	}

	double** AtWrB = new double* [8];
	for (int i = 0; i < 8; i++) {
		AtWrB[i] = new double[1];
		AtWrB[i][0] = 0;
		for (int k = 0; k < 40; k++)
			AtWrB[i][0] += AtW[i][k] * rB[k][0];
	}

	double** mTmp = new double* [8];
	for (int i = 0; i < 8; i++) {
		mTmp[i] = new double[1];
		mTmp[i][0] = 0;
		for (int k = 0; k < 8; k++)
			mTmp[i][0] += AtWA[i][k] * AtWrB[k][0];
	}

	double** rBeta = new double* [8];
	for (int i = 0; i < 8; i++) {
		rBeta[i] = new double[1];
		rBeta[i][0] = beta[i][0] - mTmp[i][0];
		B.push_back(rBeta[i][0]);
	}


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
	coeff _c = { B[4],B[5],B[6],B[7] };
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
		for (int j = 0; j < 8; ++j)
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
		km1[i] = new double[8],
		km2[i] = new double[8],
		km3[i] = new double[8],
		km4[i] = new double[8],
		stm[i] = new double[8];

		for (int j = 0; j < 8; ++j)
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
		for (int j = 0; j < 8; ++j)
			dxdp[i][j] += ((km1[i][j] + km2[i][j] * 2 + km3[i][j] * 2 + km4[i][j]) * h / 6);
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
		km1[i] = new double[8],
			km2[i] = new double[8],
			km3[i] = new double[8],
			km4[i] = new double[8],
			stm[i] = new double[8];

		for (int j = 0; j < 8; ++j)
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

	step(st, body, k2, dxdp, km2, stm, h);
	fT(st, k3, km3, stm, B);

	step(st, body, k3, dxdp, km3, stm, h * 2);
	fT(st, k4, km4, stm, B);

	for (int i = 0; i < body.size(); ++i) {
		body[i] = body[i] + ((k1[i] + k2[i] * 2 + k3[i] * 2 + k4[i]) * h / 6);
	}

	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 8; ++j)
			dxdp[i][j] += ((km1[i][j] + km2[i][j] * 2 + km3[i][j] * 2 + km4[i][j]) * h / 6);
}

bool doubleEquals(CArray arr1, CArray arr2, CArray correctArr1, CArray correctArr2) {
	double epsilon = 1006.7;
	for (int i = 0; i < 20; i++) {
		if (std::fabs(arr1[i].real() - correctArr1[i].real()) > epsilon || std::fabs(arr2[i].real() - correctArr2[i].real()) > epsilon)
			return false;
	}
	return true;
}

int main() {
	srand(time(NULL));
	double** dXdP = new double* [4];
	for (int i = 0; i < 4; i++)
	{
		dXdP[i] = new double[8];
		for (int j = 0; j < 8; ++j)
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
				for (int j = 0; j < 8; ++j)
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
	Complex correctArr1[20], correctArr2[20];
	for (int i = 0; i < 20; i++)
	{
		correctArr1[i] = correctPositions[i].first;
		correctArr2[i] = correctPositions[i].second;
	}
	Complex arr1[20], arr2[20];
	CArray data1, data2, correctData1, correctData2;
	std::vector<double> B = {150,350,0,0,5.95,6.05,0.045,0.075};
	int k = 0;
	while (true) {
		std::vector<std::pair<double, double**>> derivTmp;
		derivTmp = deriv;
		GaussNewthon(body, derivTmp, B);
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
		if (k != 0 && doubleEquals(data1, data2, correctData1, correctData2))
			break;
		double ta = 0;
		k++;
		double** dXdPT = new double* [4];
		for (int i = 0; i < 4; i++)
		{
			dXdPT[i] = new double[8];
			for (int j = 0; j < 8; ++j)
				if (i == j) dXdPT[i][j] = 1;
				else dXdPT[i][j] = 0;
		}
		std::vector<Body> bodyT;
		bodyT.push_back(Body(150, 0));
		bodyT.push_back(Body(350, 0));
		derivTmp.clear();
		positions.clear();
		while (ta <= 20) {
			RungeKuttaT(bodyT, dXdPT, h, fT,B);
			ta += h;

			if (ta - (int)ta > h) {
				positions.push_back(std::make_pair(bodyT[0].x, bodyT[1].x));
				derivTmp.push_back(std::make_pair(t, dXdPT));
			}
		}
	}
	std::cout << B[4] << ' ' << B[5] << ' ' << B[6] << ' ' << B[7] << '\n';
	return 0;
}
