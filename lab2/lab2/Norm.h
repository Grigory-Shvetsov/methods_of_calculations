#pragma once
using namespace std;
#include<vector>
#include<fstream>

template<typename T>
using matrix = vector<vector<T>>;



// норма p=2 квадратной матрицы
template <typename T>
T normMatr(const matrix<T>& A)
{
	size_t n = A.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < n; ++j)
			res += A[i][j] * A[i][j];

	res = sqrt(res);
	return res;
}

// норма p=1 квадратной матрицы
template <typename T>
T norm1Matr(const matrix<T>& A)
{
	size_t n = A.size();
	T max = 0;
	for (size_t j = 0; j < n; ++j)
	{
		T sum = 0;
		for (size_t i = 0; i < n; ++i)
			sum += abs(A[i][j]);
		if (sum > max)
			max = sum;
	}
	return max;
}

// норма p=Inf квадратной матрицы
template <typename T>
T normInfMatr(const matrix<T>& A)
{
	size_t n = A.size();
	T max = 0;
	for (size_t i = 0; i < n; ++i)
	{
		T sum = 0;
		for (size_t j = 0; j < n; ++j)
			sum += abs(A[i][j]);
		if (sum > max)
			max = sum;
	}

	return max;
}


// норма p=2 вектора
template <typename T>
T norm(const vector<T>& A)
{
	size_t n = A.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		res += A[i] * A[i];

	res = sqrt(res);
	return res;
}

// норма p=1 вектора
template <typename T>
T norm1(const vector<T>& A)
{
	size_t n = A.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		res += abs(A[i]);

	return res;
}

// норма p=Inf вектора
template <typename T>
T normInf(const vector<T>& A)
{
	size_t n = A.size();
	T max = 0;
	for (size_t i = 0; i < n; ++i)
		if (abs(A[i]) > max)
			max = abs(A[i]);

	return max;
}



//¬ывод матрицы ј
template<typename T>
void print_vec(const matrix<T>& vec_A)
{
	for (auto& i : vec_A) {
		for (auto j : i)
			cout << j << "\t";
		cout << endl;
	}
}


//¬ывод расширенной матрицы
template<typename T>
void print_AB(const matrix<T>& vec_A, const vector<T>& vec_B)
{
	cout << endl;
	for (size_t i = 0; i < vec_B.size(); i++) {
		for (size_t j = 0; j < vec_B.size(); j++)
			cout << fixed << setprecision(7) << vec_A[i][j] << "\t";
		cout << fixed << setprecision(7) << vec_B[i] << endl;
	}
}



//¬ывод вектора
template<typename T>
void print_vec(const vector<T>& vec)
{
	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << fixed << setprecision(8) << vec[i] << "\t";
	cout << endl;
}

//считывание из файла
template<typename T>
size_t out_of_file(matrix<T>& vec_A, vector<T>& vec_B, string str)
{
	size_t N = 0;
	ifstream fout(str);
	fout >> N;
	T x_A, x_B;
	for (size_t i = 0; i < N; ++i) {
		vector<T> v_A;
		for (size_t j = 0; j < N; ++j) {
			fout >> x_A;
			v_A.push_back(x_A);
		}
		vec_A.push_back(v_A);
		fout >> x_B;
		vec_B.push_back(x_B);
	}
	fout.close();
	return N;
}

// умножение матрицы и столбца
template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n); // оптимизируем пам€ть
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

// умножение матриц
template<typename T>
matrix<T> mult(const matrix<T>& A1, const matrix<T>& A2) {
	size_t l = A1.size();
	size_t m = A2.size();
	size_t n = A2[0].size();

	vector<vector<T>> res;
	res.reserve(l);
	for (size_t i = 0; i < l; ++i) {
		vector<T> row;
		row.reserve(n); // оптимизируем пам€ть
		for (size_t j = 0; j < n; ++j) {
			T c = 0;
			for (size_t k = 0; k < m; ++k)
				c += A1[i][k] * A2[k][j];
			row.push_back(c);
		}
		res.push_back(row);
	}
	return res;
}


// —оздает единичную матрицу размером NxN
template<typename T>
matrix<T> eig(size_t n) {
	vector<vector<T>> e;
	e.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row(n, 0);
		row[i] = 1;
		e.push_back(row);
	}
	return e;
}

// транспонирование квадратной матрицы (измен€ет входную матрицу на транспорированую)
template<typename T>
void transpose(matrix<T>& A) {
	size_t n = A.size();
	// тут была ошибка
	for (size_t i = 0; i < n - 1; ++i)
		for (size_t j = i + 1; j < n; ++j)
			swap(A[i][j], A[j][i]);
}

//создает нулевую матрицу N*N
template<typename T>
matrix<T> zero(size_t N)
{
	matrix<T> res;
	for (size_t i = 0; i < N; i++) {
		vector<T> row(N, 0);
		res.push_back(row);
	}
	return res;
}

//создает трехдиагональную матрицу
template<typename T>
pair<matrix<T>, vector<T>> test3(size_t N)
{
	matrix<T> A = zero<T>(200 + N);
	vector<T> B;
	B.reserve(200 + N);
	A[0][0] = 4;
	A[0][1] = 1;
	for (size_t i = 1; i < 200 + N - 1; i++) {
		A[i][i] = 4;
		A[i][i + 1] = 1;
		A[i][i - 1] = 1;
	}
	A[200 + N - 1][200 + N - 1] = 4;
	A[200 + N - 1][200 + N - 2] = 1;


	B.push_back(6);
	for (size_t i = 1; i < 200 + N - 1; i++) {
		// прибавл€ем единичку к i из-за индексации с нул€
		B.push_back(10 - 2 * ((i + 1) % 2));
	}
	B.push_back(9 - 3 * ((200 + N) % 2));

	return make_pair(A, B);
}

// создает трехдиагональную матрицу дл€ 3 теста
// ¬озвращает:
// 1. нижнюю диагональ
// 2. главную диагональ
// 3. верхнюю диагональ
// 4. вектор правой части
template<typename T>
tuple<vector<T>, vector<T>, vector<T>, vector<T>> test3d(size_t N) {
	vector<T> a(200 + N - 1, 1);
	vector<T> b(200 + N, 4);
	vector<T> c(200 + N - 1, 1);

	vector<T> B;
	B.reserve(200 + N);
	B.push_back(6);
	for (size_t i = 1; i < 200 + N - 1; i++) {
		// прибавл€ем единичку к i из-за индексации с нул€
		B.push_back(10 - 2 * ((i + 1) % 2));
	}
	B.push_back(9 - 3 * ((200 + N) % 2));

	return make_tuple(a, b, c, B);
}

////
// перегруженные операторы
////

// умножение матрицы на вектор
template<typename T>
vector<T> operator*(const matrix<T>& A1, const vector<T>& A2) {
	return mult(A1, A2);
}

// умножение двух матриц
template<typename T>
matrix<T> operator*(const matrix<T>& A1, const matrix<T>& A2) {
	return mult(A1, A2);
}

// умножение вектора на скал€р без выделени€ пам€ти
template<typename T>
vector<T>& operator*(vector<T>&& v, T b) {
	for (auto& el : v)
		el *= b;
	return v;
}

// умножение вектора на скал€р
template<typename T>
vector<T> operator*(const vector<T>& v, T b) {
	vector<T> res(v);
	return move(res) * b;
}

// умножение матрицы на скал€р без выделени€ пам€ти
template<typename T>
matrix<T>& operator*(matrix<T>&& A, T b) {
	for (auto& row : A)
		for (auto& el : row)
			el *= b;
	return A;
}

// умножение матрицы на скал€р
template<typename T>
matrix<T> operator*(const matrix<T>& A, T b) {
	return matrix<T>(A) * b;
}

// сложение векторов без выделени€ пам€ти
template<typename T>
vector<T>& operator+(vector<T>&& A1, const vector<T>& A2) {
	size_t n = A1.size();
	for (size_t i = 0; i < n; ++i)
		A1[i] += A2[i];
	return A1;
}

// сложение векторов
template<typename T>
vector<T> operator+(const vector<T>& A1, const vector<T>& A2) {
	vector<T> res(A1);
	return move(res) + A2;
}

// сложение матриц без выделени€ пам€ти
template<typename T>
matrix<T>& operator+(matrix<T>&& A1, const matrix<T>& A2) {
	size_t n = A1.size();
	size_t m = A1[0].size();
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			A1[i][j] += A2[i][j];
	return A1;
}

// сложение матриц
template<typename T>
matrix<T> operator+(const matrix<T>& A1, const matrix<T>& A2) {
	matrix<T> res(A1);
	return move(res) + A2;
}

// вычитание векторов без выделени€ пам€ти
template<typename T>
vector<T>& operator-(vector<T>&& A1, const vector<T>& A2) {
	size_t n = A1.size();
	for (size_t i = 0; i < n; ++i)
		A1[i] -= A2[i];
	return A1;
}

// вычитание векторов
template<typename T>
vector<T> operator-(const vector<T>& A1, const vector<T>& A2) {
	vector<T> res(A1);
	return move(res) - A2;
}

// вычитание матриц без выделени€ пам€ти
template<typename T>
matrix<T>& operator-(matrix<T>&& A1, const matrix<T>& A2) {
	size_t n = A1.size();
	size_t m = A1[0].size();
	for (size_t i = 0; i < n; ++i)
		for (size_t j = 0; j < m; ++j)
			A1[i][j] -= A2[i][j];
	return A1;
}

// сложение матриц
template<typename T>
matrix<T> operator-(const matrix<T>& A1, const matrix<T>& A2) {
	matrix<T> res(A1);
	return move(res) - A2;
}


//ѕрава€ часть
template<typename T>
vector<T> right_part(const vector<vector<T>>& vec_A, const vector<T>& vec_X)
{
	T b1;
	vector<T> vec_B1;
	for (size_t i = 0; i < vec_X.size(); i++) {
		b1 = 0;
		for (size_t j = 0; j < vec_X.size(); j++)
			b1 += vec_A[i][j] * vec_X[j];
		vec_B1.push_back(b1);
	}
	return vec_B1;
}



//нев€зка
template<typename T>
T nevazka(const matrix<T>& vec_A_source, const vector<T>& vec_B_source, const vector<T>& vec_X)
{
	vector<T> vec_AX = right_part(vec_A_source, vec_X);
	vector<T> vec;
	T bb1 = 0;
	for (size_t i = 0; i < vec_AX.size(); i++) {
		bb1 = vec_AX[i] - vec_B_source[i];
		vec.push_back(bb1);
	}
	return normInf(vec);
}

