#pragma once
#include <iostream>
#include <vector>

using namespace std;
template<typename T>
using matrix = vector<vector<T>>;


//Вывод вектора
template<typename T>
void print_vec(const vector<T>& vec)
{

	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << vec[i] << "\t";
	cout << endl << endl;
}

//Вывод вектора пар локализации
template<typename T>
void print_vec(const vector<pair<T, T>>& vec)
{
	for (size_t i = 0; i < vec.size(); i++)
		cout << vec[i].first << "\t" << vec[i].second << endl;
	cout << endl;
}



// умножение матрицы и столбца
template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n); // оптимизируем память
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

// норма (евклидова) вектора
template<typename T>
T norm(const vector<T>& x) {
	T res = 0;
	size_t n = x.size();
	for (size_t i = 0; i < n; ++i)
		res += x[i] * x[i];
	return sqrt(res);
}


// Обращение матрицы 2x2 (просто ввел формулы из Wolfram Mathematica)
template<typename T>
matrix<T> inv(const matrix<T>& m) {
	T det = m[0][0] * m[1][1] - m[1][0] * m[0][1];
	return matrix<T> {
		{  m[1][1] / det, -m[0][1] / det },
		{ -m[1][0] / det,  m[0][0] / det }
	};
}



 // ошибка
 template<typename T>
 void degreeX(vector<T> X, T x)
 {
	 for (size_t i = 0; i < X.size()-2; i++) {
		 cout << "Порядок p = " << log(abs((X[i+2] - x) / (X[i+1] - x))) / log(abs((X[i+1] - x) / (X[i] - x))) << endl;
	 }
 }

// Численная производная функции
template<typename T, typename F>
auto numDeriv(F& fun, T eps) {
	return [&fun, eps](T x) {
		return (fun(x + eps) - fun(x)) / eps;
	};
}

// Вывод отрезков локализации
template<typename T>
ostream& operator<<(ostream& out, vector<pair<T, T>> intervals) {
	if (intervals.size() > 0) {
		pair<T, T> inter = intervals[0];
		out << '[' << inter.first << ", " << inter.second << "]";
	}
	for (size_t i = 1; i < intervals.size(); ++i) {
		pair<T, T> inter = intervals[i];
		out << "; [" << inter.first << ", " << inter.second << "]";
	}
	return out;
}