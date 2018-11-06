#ifndef _complex_cpp
#define _complex_cpp
#include "vector.h"
#include<iostream> 
#define M_PI            3.14159265358979323846
#define M_PI_2          1.57079632679489661923

//конст+дестр
vector::vector() : N(0){}
vector::vector(unsigned n) : N(n) {
	A = new double[N];
	for (unsigned i = 0; i < N; ++i)
		A[i] = 0;
}
vector::vector(unsigned n, double a) : N(n) {
	A = new double[N];
	for (unsigned i = 0; i < N; ++i)
		A[i] = a;
}
vector::vector(const vector&z) : N(z.N)
{
	A = new double[N];
	for (unsigned i = 0; i < N; ++i)
		A[i] = z.A[i];
} 
vector::~vector() { delete[]A; }




const vector& vector :: operator=(const vector&z)
{
	if (N <= z.N && N!=0) 
		memcpy(A,z.A, sizeof(*z.A) * N);
	else if (N == 0)
	{
		if (z.N > 0)
		{
			A = new double[N=z.N];
			memcpy(A, z.A, sizeof(*z.A)*N);
		}
	}
	else if (N > z.N)
	{
		memcpy(A, z.A, sizeof(*z.A) * z.N);
		for (unsigned i = z.N; i < N; ++i)
			A[i] = 0;
	}
	return *this;
}
double & vector:: operator[](unsigned num)const 
{
	if (num >= N)
		std::cout << "Номер индекса больше количества элементов." << std::endl;
	else
		return A[num];
}
const vector vector :: operator+() { return *this; }
const vector vector :: operator-()
{
	vector z2(N);
	if (N > 0) 
	{
		for (unsigned i = 0; i < N; ++i)
			z2.A[i] = -A[i];
		
		return z2;
	}
	else 
		return *this;
}
const vector& vector :: operator++()const
{
	if (N > 0)
	{
		for (unsigned i = 0; i < N; ++i)
			A[i] += 1;
	}
	return *this;
}
const vector& vector :: operator--()const
{
	if (N > 0)
	{
		for (unsigned i = 0; i < N; ++i)
			A[i] -= 1;

	}
	return *this;
}
const vector vector :: operator++(int)const
{
	vector a = *this;
	if (N > 0)
	{
		for (unsigned i = 0; i < N; ++i)
			A[i] += 1;
	}
	return a;
}
const vector vector :: operator--(int)const
{
	vector a = *this;
	if (N > 0)
	{
		for (unsigned i = 0; i < N; ++i)
			A[i] -= 1;
	}
	return a;
}

const bool vector :: operator ==(const vector &z)
{
	if (N == z.N)
	{
		if (memcmp(A, z.A, sizeof(*z.A) * N) == 0)
		{
			return true;
		}
		else
			return false;
	}
	else 
		return false;
}
const bool vector :: operator !=(const vector &z)
{
	if (N == z.N)
	{
		if (memcmp(A, z.A, sizeof(*z.A) * N) == 0)
		{
			return false;
		}
		else
			return true;
	}
	else
		return true;
}



const vector& vector :: operator+= (const vector &z)
{
	if (z.N >= N)
	{
		for(unsigned i = 0; i < N; ++i)
			A[i] += z.A[i];
	}
	else 
		if (z.N <= N) 
		{
			for (unsigned i = 0; i <z.N; ++i)
				A[i] += z.A[i];
		}
	return *this;
}
const vector& vector :: operator-= (const vector &z)
{
	if (z.N >= N)
	{
		for (unsigned i = 0; i < N; ++i)
			A[i] -= z.A[i];
	}
	else
		if (z.N < N)
		{
			for (unsigned i = 0; i < z.N; ++i)
				A[i] -= z.A[i];
		}
	return *this;
}

const vector& vector :: operator*= (const double &z)
{
	for (unsigned i = 0; i < N; ++i)
		A[i] *= z;
	return *this;
}
const vector& vector :: operator/= (const double &z)
{
	for (unsigned i = 0; i < N; ++i)
		A[i] /= z;
	return *this;
}


const vector vector :: operator+ (const vector &z) const {
	vector z2;
	if (N <= z.N)
	{
		z2.A = new double[z2.N = z.N];
		for (unsigned i = 0; i < N; ++i)
			z2.A[i] = A[i]+ z.A[i];
		for (unsigned i = N; i < z.N; ++i)
			z2.A[i] = z.A[i];

	}
	else 
	{
		z2.A = new double[z2.N = N];
		for (unsigned i = 0; i < z.N; ++i)
			z2.A[i] = A[i] + z.A[i];
		for (unsigned i = z.N; i < N; ++i)
			z2.A[i] = A[i];
	}
	return z2;
}
const vector vector :: operator- (const vector &z) const 
{
	vector z2;
	if (N <= z.N)
	{
		z2.A = new double[z2.N = z.N];
		for (unsigned i = 0; i < N; ++i)
			z2.A[i] = A[i] - z.A[i];
		for (unsigned i = N; i < z.N; ++i)
			z2.A[i] = z.A[i];

	}
	else
	{
		z2.A = new double[z2.N = N];
		for (unsigned i = 0; i < z.N; ++i)
			z2.A[i] = A[i] - z.A[i];
		for (unsigned i = z.N; i < N; ++i)
			z2.A[i] = A[i];
	}
	return z2;
}


const double vector :: operator* (const vector &z) const 
{ 
	double a = 0;
	if (N == z.N)
	{
		for (unsigned i = 0; i < N; ++i)
			a += A[i] * z.A[i];
	}
	else std::cout << "Размерность векторов не совпадает"<<std::endl;
	return a;
}

const vector vector :: operator* (const double &z) const
{
	vector a(N);
		for (unsigned i = 0; i < N; ++i)
			a.A[i] = A[i] * z;
	return a;
}
const vector vector :: operator/ (const double &z) const
{
	vector a(N);
	for (unsigned i = 0; i < N; ++i)
		a.A[i] = A[i] / z;
	return a;
}


vector operator *(double c, const vector&z)
{
	vector a(z.N);
	for (unsigned i = 0; i < z.N; ++i)
		a.A[i] = c*z.A[i];
	return a;
}

//длинна ветора
const double vector :: len() const
{
	double a=0;
	for (unsigned i = 0; i < N; ++i)
		a += A[i] * A[i];
	return sqrt(a);
}


//<< и >>
std::ostream & operator<<(std::ostream &c, const vector &z) {
	c << "Размерность---" << z.N << std::endl << "(";
	for (unsigned i = 0; i < z.N; ++i)
	{

		c << z.A[i];
		if (i != z.N - 1) c << ", ";
	}
	c << ")";
	return c;

}
std::istream & operator>>(std::istream &c, vector &z) {
	double buf_size = 4;
	double value;
	double *buf = new double[buf_size];
	unsigned q = 0;
	for (;std::cin >> value; ++q)
	{
		if (q >= buf_size)
		{
			double nb =q+1;
			double *n_buf = new double[nb];
			for (int u = 0; u < q; ++u)
				n_buf[u] = buf[u];
			delete[] buf;
			buf = n_buf;
			buf_size = nb;
		}
		buf[q] = value;
		std::cout << buf[q]<<std::endl;
	}
	z.N = q;
	z.A = new double[z.N];
	for (unsigned i = 0; i < z.N; ++i)
		z.A[i]= buf[i];
	return c;

}

unsigned vector :: getN() { return N;}
double vector::getAi(unsigned i)
{
	if (i >= N)
	{
		std::cout << "Нет такого индекса" << std::endl;
		return 0;
	}
	else return A[i];
}

void vector::setN(unsigned i) 
{ 
	N = i;
	A = new double[N];
}
void vector::setAi(unsigned i, double c)
{

	if (i >= N)
	{
		std::cout << "Нет такого индекса" << std::endl;
	}
	else A[i]=c;
}






/*double vector::setAi(unsigned i)


{
	return A[i];
}

/*
const complex complex :: operator/ (const complex &z) const { return complex::div(z); }
const complex complex :: operator+ (const double &z) const
{
	complex z1;
	z1.re = re + z;
	z1.im = im;
	return z1;
}
const complex complex :: operator- (const double &z) const
{
	complex z1;
	z1.re = re - z;
	z1.im = im;
	return z1;
}
const complex complex :: operator* (const double &z) const
{
	complex z1;
	z1.re = re * z;
	z1.im = im * z;
	return z1;
}
const complex complex :: operator/ (const double &z) const
{
	complex z1;
	z1.re = re / z;
	z1.im = im / z;
	return z1;
}


const complex complex :: operator+= (const complex &z)
{
	re=re + z.re;
	im=im + z.im;
	return *this;
}
const complex complex :: operator-= (const complex &z)
{
	re = re - z.re;
	im = im - z.im;
	return *this;
}
const complex complex :: operator*= (const complex &z)
{
	double r=re;
	re = re * z.re - im * z.im;
	im = r * z.im + im * z.re;
	return *this;
}
const complex complex :: operator/= (const complex &z)
{
	double r=re;
	re = (re*z.re + im * z.im) / (z.re*z.re + z.im*z.im);
	im = (z.re*im - r * z.im) / (z.re*z.re + z.im*z.im);
	return *this;
}



const complex complex :: operator+= (const double &z)
{
	re += z;
	return *this;
}
const complex complex :: operator-= (const double &z)
{
	re -= z;
	return *this;
}
const complex complex :: operator*= (const double &z)
{
	re *= z;
	im *= z;
	return *this;
}
const complex complex :: operator/= (const double &z)
{
	re /= z;
	im /= z;
	return *this;
}


const complex complex::step(double st) const
{
	double r = sqrt(re*re + im * im), u = fi();
	complex ot;
	ot.setRe(pow(r, st)*(cos(u*st)));
	ot.setIm(pow(r, st)*(sin(u*st)));
	return ot;
}



complex::complex(double Re, double Im) : re(Re), im(Im) { }
complex::complex() : re(0), im(0) { }
std::ostream & operator<<(std::ostream &c, const complex &z) {
	if (z.im > 0) { c << z.re << "+" << z.im << "i"; }
	else
	{
		if (z.im < 0) { c << z.re << z.im << "i"; }
		else c << z.re;
	}
	return c;

}
std::istream & operator>>(std::istream &c, complex &z) {
	c >> z.re >> z.im;
	return c;

}
const complex operator+(double re, const complex &z) {
	return z + re;
}
const complex operator-(double re, const complex &z) {
	complex z1;
	z1.re = re - z.re;
	z1.im = -z.im;
	return z1;;
}
const complex operator*(double re, const complex &z) {
	return z * re;
}
const complex operator/(double re, const complex &z) {
	complex z1;
	z1.re = (re * z.re) / (z.re*z.re + z.im*z.im);
	z1.im = ((-z.im)*re) / (z.re*z.re + z.im*z.im);
	return z1;;
}
const double complex::fi() const
{
	if (re > 0 && im >= 0)return std::atan(im / re);
	else if (re < 0 && im >= 0) return M_PI - atan(fabs(im / re));
	else if (re < 0 && im < 0) return M_PI + atan(fabs(im / re));
	else if (re > 0 && im < 0) return M_PI + M_PI - atan(fabs(im / re));
	else if (re == 0 && im > 0) return M_PI_2;
	else if (re == 0 && im < 0) return 3 * M_PI_2;
}
void complex::trig()
{
	std::cout << sqrt(re*re + im * im) << "(cos(" << fi() << ")+sin(" << fi() << ")";
}
void complex::setRe(double Re) { re = Re; }
void complex::setIm(double Im) { im = Im; }
double complex::getRe() { return re; }
double complex::getIm() { return im; }
*/
#endif 