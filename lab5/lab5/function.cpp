#pragma once
#include <iostream>
#include <vector>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;

// TODO добавить производные примеров




//0 <= x <= 1
template<typename T>
T test0(T x)
{
	return (x - 1) * (x - 1) * (x - 3) * (x - 3);
}

template<typename T>
T test0Df(T x)
{
	return 4 * (-6 + 11 * x - 6 * x * x + x * x * x);
}




//0 <= x <= 1
template<typename T>
T test1(T x)
{
	return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
}

template<typename T>
T test1Df(T x)
{
	return 0.121495 - 1.5119 * x + 5.9535 * x * x - 9.28 * x * x * x + 5. * x * x * x * x;
}


//-1 <= x <= 10
template<typename T>
T test2(T x)
{
	return sqrt(x + 1) - 1;
}

template<typename T>
T test2Df(T x)
{
	return 1. / (2. * sqrt(1. + x));
}

//0 <= x <= 1
template<typename T>
T test3(T x)
{
	return 35 * x * x * x - 67 * x * x - 3 * x + 3;
}

template<typename T>
T test3Df(T x)
{
	return -3. - 134. * x + 105. * x * x;
}

// |x|<=10, |y|<=10
template<typename T>
pair<T, T> test4(T x, T y)
{
	return make_pair(x * x - y * y - 15, x * y + 4);
}

template<typename T>
T test4f1(T x, T y)
{
	return x * x - y * y - 15.;
}

template<typename T>
T test4f2(T x, T y)
{
	return x * y + 4.;
}

template<typename T>
matrix<T> test4Jacobi(T x, T y)
{
	return matrix<T>{
		{2. * x, -2. * y},
		{ y, x }
	};
}

// |x|<=10, |y|<=10
template<typename T>
pair<T, T> test5(T x, T y)
{
	return make_pair(x * x + y * y + x + y - 8, x * x + y * y + x * y - 7);
}

template<typename T>
T test5f1(T x, T y)
{
	return x * x + y * y + x + y - 8.;
}

template<typename T>
T test5f2(T x, T y)
{
	return x * x + y * y + x * y - 7.;
}

template<typename T>
matrix<T> test5Jacobi(T x, T y)
{
	return matrix<T>{
		{1. + 2. * x, 1. + 2. * y},
		{ 2. * x + y, x + 2. * y }
	};
}


template<typename T>
T test6f1(T x, T y)
{
	return x * x + y * y - 1.;
}

template<typename T>
T test6f2(T x, T y)
{
	return x * x * x - y;
}

template<typename T>
matrix<T> test6Jacobi(T x, T y)
{
	return matrix<T>{
		{2. * x, 2. * y},
		{ 3. * x * x, -1. }
	};
}


// Тест для отчёта
// -1 <= x <= 1
// Один корень приблизительно равный 0.516
template<typename T>
T test7(T x) {
	return x * x - 4 * sin(x) + log(x + 5.);
}
template<typename T>
T test7Df(T x) {
	return 2. * x - 4 * cos(x) + 1. / (x + 5.);
}



