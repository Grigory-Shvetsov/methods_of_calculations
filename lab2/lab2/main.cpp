#include<iostream>
#include<iomanip>
#include<string>
#include<algorithm>
#include"Norm.h"


vector<double>Test1{ 5,-7,12,4 };
vector<double>Test2{ 10,-10,12,4 };


template<typename T>
bool stop(const vector<T>& Xk, const vector<T>& Xk1, T coef, T eps)
{
	return normInf(Xk1 - Xk) > coef;
}


//меняет строки матрицы А и В чтобы норма стала меньше единицы
template<typename T>
void myswap(matrix<T>& A, vector<T>& B, T tau)
{
	T normC = normInfMatr(eig<T>(B.size()) - A * tau);
	if (normC <= 1.0) {}
	else
	{
		cout << "Норма перед преобразованиями ||C|| = " << normC << endl;
		for (size_t i = 0; i < B.size() - 1; i++)
			for (size_t j = i + 1; j < B.size(); j++)
			{
				cout << "Меняются строки" << endl;
				swap(A[i], A[j]);
				swap(B[i], B[j]);
				print_AB(A, B);
				cout << "||C|| = " << normInfMatr(eig<T>(B.size()) - A * tau) << endl;
				if (normInfMatr(eig<T>(B.size()) - A * tau) <= 1.0)
					break;
				if (normC <= normInfMatr(eig<T>(B.size()) - A * tau)) {
					swap(A[i], A[j]);
					swap(B[i], B[j]);
					print_AB(A, B);
					cout << "________________________________________________________________" << endl;
				}
				normC = normInfMatr(eig<T>(B.size()) - A * tau);
			}
	}
}


//меняет строки матрицы А и В чтобы норма стала меньше единицы
template<typename T>
void myminus(matrix<T>& A, vector<T>& B, T tau)
{
	T normC = normInfMatr(eig<T>(B.size()) - A * tau);
	if (normC <= 1.0) {}
	else
	{
		cout << "Норма перед преобразованиями ||C|| = " << normC << endl;
		for (size_t i = 0; i < B.size() - 1; i++)
		{
			A[i] = operator*(A[i], -1.0);
			B[i] = B[i] * (-1);
			for (size_t j = i + 1; j < B.size(); j++)
			{
				T normC1 = normInfMatr(eig<T>(B.size()) - A * tau);
				A[j] = operator*(A[j], -1.0);
				B[j] = B[j] * (-1);
				T normC2 = normInfMatr(eig<T>(B.size()) - A * tau);

				if (normC1 < normC2) {
					A[j] = operator*(A[j], -1.0);
					B[j] = B[j] * (-1);
					normC = normC1;
				}
				else
					normC = normC2;
				print_AB(A, B);
				cout << "||C|| = " << normC << endl;


				if (normInfMatr(eig<T>(B.size()) - A * tau) < 1)
					break;
				normC = normInfMatr(eig<T>(B.size()) - A * tau);

			}

		}
	}
}

// Метод простой итерации
// X0 - начальная точка
// tau - параметр метода
// eps - погрешность
// По заданию (пункт 7) tau нужно выбирать экспериментально


//оценка для числа итераций
//r0 - метрика(начального приближения, следующего элемента)
int k_est(double normC, double r0, double eps)
{
	return ceil(log((1 - normC) * eps / r0) / log(normC));
}

template<typename T>
vector<T> simpleIterationMethod(
	matrix<T>& A,
	vector<T>& B,
	const vector<T>& X0,
	T tau,
	T eps
) {
	int k = 0;
	size_t n = A.size();
	//myswap(A, B, tau); // вдруг получится изменить норму перестановкой строк
	//myminus(A, B, tau); // вдруг получится изменить норму домножением на -1
	matrix<T> C = eig<T>(n) - A * tau;
	T normC = normInfMatr(C);
	cout << endl << "||C|| = " << normC << endl;
	/*cout << endl << "Матрица С:";
	print_vec(C);*/

	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;

	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}




	vector<T> tauB = B * tau;
	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1;

	auto kEst = k_est(normC, normInf(C * Xk + tauB - Xk), eps);
	do
	{
		Xk1 = C * Xk + tauB;
		Xk_ = Xk;
		Xk = Xk1;
		++k;
	} while (stop(Xk_, Xk1, coef, eps));
	cout << "k = " << k << " <= " << kEst << endl;

	vector<T> Xk_est(X0);
	vector<T> Xk1_est;
	for (size_t i = 0; i < kEst; i++) {
		Xk1_est = C * Xk_est + tauB;
		Xk_est = Xk1_est;
	}
	print_vec(Xk1_est);
	cout << "Норма ошибки после k_est = " << normInf(Xk1_est - Test2) << endl;
	cout << "|| Ошибка || = " << normInf(Xk1 - Test2) << endl;

	return Xk1;
}

// Метод Якоби
template<typename T>
vector<T> JacobiMethod(
	const matrix<T>& A,
	const vector<T>& B,
	const vector<T>& X0,
	T eps
) {
	int k = 0;
	size_t n = A.size();
	matrix<T> C;
	C.reserve(n);
	vector<T> B1;
	B1.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		T Aii = A[i][i];
		// Заполняем C
		vector<T> rowC;
		rowC.reserve(n);
		for (size_t j = 0; j < i; ++j)
			rowC.push_back(-A[i][j] / Aii);

		rowC.push_back(0);

		for (size_t j = i + 1; j < n; ++j)
			rowC.push_back(-A[i][j] / Aii);

		C.push_back(rowC);

		// Заполняем B1 ("y" в методичке)
		B1.push_back(B[i] / Aii);
	}
	T normC = normInfMatr(C);
	cout << "||C|| = " << fixed << setprecision(8) << normC << endl;
	/*cout << endl << "Матрица С:";
	print_vec(C);*/
	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;

	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}

	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1;
	// норма невязки
	auto kEst = k_est(normC, normInf(C * Xk + B1 - Xk), eps);
	do
	{
		Xk1 = C * Xk + B1;
		Xk_ = Xk;
		Xk = Xk1;
		++k;
	} while (nevazka(A, B, Xk1) > eps);
	cout << "k = " << k << " <= " << kEst << endl;

	vector<T> Xk_est(X0);
	vector<T> Xk1_est;
	for (size_t i = 0; i < kEst; i++) {
		Xk1_est = C * Xk_est + B1;
		Xk_est = Xk1_est;
	}
	print_vec(Xk1_est);
	cout << "Норма ошибки после k_est = " << normInf(Xk1_est - Test2) << endl;
	cout << "|| Ошибка || = " << normInf(Xk1 - Test2) << endl;
	return Xk1;

}

// Обратная к нижнетругольной матрице (заполняет строки не до конца,
// результирующая матрица будет нижнетреугольной).
// Использутеся для вычисления матрицы C в методах Зейделя и релаксации
template<typename T>
matrix<T> inverseL(const matrix<T>& L) {
	size_t n = L.size();
	matrix<T> res;
	res.reserve(n);

	vector<T> row1{ 1 / L[0][0] };
	res.push_back(row1);
	for (size_t i = 1; i < n; ++i) {
		vector<T> row;
		row.reserve(i + 1);

		T resii = 1 / L[i][i];
		for (size_t j = 0; j < i; ++j) {
			T sum = 0;
			for (size_t k = j; k < i; ++k)
				sum += res[k][j] * L[i][k];
			row.push_back(-resii * sum);
		}
		row.push_back(resii);
		res.push_back(row);
	}
	return res;
}

// Произведение нижнетреугольной и верхнетреугольной матриц.
// Метод предполагает, что элементы под диагональю верхнетреугольной матрицы не
// заполнены (вместо них идут сразу диагольный элемент и т.д.)
// Использутеся для вычисления матрицы C в методах Зейделя и релаксации.
template<typename T>
matrix<T> multLU(const matrix<T>& L, const matrix<T>& U) {
	size_t n = L.size();
	matrix<T> res;
	res.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row;
		row.reserve(n);

		for (size_t j = 0; j < n; ++j) {
			size_t m = min(i, j);

			T sum = 0;
			for (size_t k = 0; k <= m; ++k)
				//sum += L[i][k] * U[k][j];
				// учитываем отсутствующие элементы верхнетреугольной матрицы
				sum += L[i][k] * U[k][j - k];
			row.push_back(sum);
		}
		res.push_back(row);
	}
	return res;
}

// Возвращает матрицу С для методов Зейделя (omega = 1) и релаксации
template<typename T>
matrix<T> relaxC(const matrix<T>& A, T omega) {
	size_t n = A.size();
	// C = (D + wL)^-1 ((1-w)D - wU)
	// D+wL
	// xk1 = Cxk + y
	matrix<T> DwL;
	DwL.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row;
		row.reserve(i + 1);

		for (size_t j = 0; j < i; ++j)
			row.push_back(omega * A[i][j]);
		row.push_back(A[i][i]);
		DwL.push_back(row);
	}
	matrix<T> DwLinv = inverseL(DwL);

	// (1-w)D - wU
	matrix<T> wDwU;
	wDwU.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row;
		row.reserve(i + 1);

		// (1-w)D
		row.push_back((1 - omega) * A[i][i]);

		// -wU
		for (size_t j = i + 1; j < n; ++j)
			row.push_back(-omega * A[i][j]);

		wDwU.push_back(row);
	}

	return multLU(DwLinv, wDwU);
}

// Метод Зейделя
template<typename T>
vector<T> SeidelMethod(
	const matrix<T>& A,
	const vector<T>& B,
	const vector<T>& X0,
	T eps
) {
	int k = 0;
	size_t n = A.size();
	matrix<T> C = relaxC(A, T(1));
	T normC = normInfMatr(C);
	cout << "||C|| = " << normC << endl;
	/*cout << endl << "Матрица С:";
	print_vec(C);*/
	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1(n, 0);

	// коэффициент для критерия выхода
	//T coef = eps;
	T coef = (1 - normC) / normC * eps;

	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}
	int kEst = 0;
	do
	{
		for (size_t i = 0; i < n; ++i) {
			// сумма со слагаемыми X^(k+1)
			T sumXk1 = 0;
			// проверка, чтобы не было переполнения при вычислении i-1
			if (i != 0)
				for (size_t j = 0; j < i; ++j)
					sumXk1 += A[i][j] * Xk1[j];

			// сумма со слагаемыми X^k (от прошлой итерации)
			T sumXk = 0;
			for (size_t j = i + 1; j < n; ++j)
				sumXk += A[i][j] * Xk[j];

			Xk1[i] = (B[i] - sumXk1 - sumXk) / A[i][i];
		}
		Xk_ = Xk;
		Xk = Xk1;
		++k;
		if (k == 1)
			kEst = k_est(normC, normInf(Xk1 - Xk_), eps);

	} while (nevazka(A, B, Xk1) > eps);

	cout << "k = " << k << " <= " << kEst << endl;
	cout << "|| Ошибка || = " << normInf(Xk1 - Test2) << endl;

	return Xk1;
}

// Релаксационный метод
// omega - параметр релаксации; omega > 0. Для симметричных матриц метод
// сходится при 0 < omega < 2
template<typename T>
vector<T> RelaxMethod(
	const matrix<T>& A,
	const vector<T>& B,
	const vector<T>& X0,
	T omega,
	T eps
) {
	int k = 0;
	size_t n = A.size();
	matrix<T> C = relaxC(A, omega);
	T normC = normInfMatr(C);
	cout << "||C|| = " << normC << endl;
	/*cout << endl << "Матрица С:";
	print_vec(C);*/
	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1(n, 0);
	// норма невязки


	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;

	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}
	int kEst = 0;
	do
	{
		for (size_t i = 0; i < n; ++i) {
			// сумма со слагаемыми X^(k+1)
			T sumXk1 = 0;
			// проверка, чтобы не было переполнения при вычислении i-1
			if (i != 0)
				for (size_t j = 0; j < i; ++j)
					sumXk1 += A[i][j] * Xk1[j];

			// сумма со слагаемыми X^k (от прошлой итерации)
			T sumXk = 0;
			for (size_t j = i + 1; j < n; ++j)
				sumXk += A[i][j] * Xk[j];

			Xk1[i] = omega * (B[i] - sumXk1 - sumXk) / A[i][i] + (1 - omega) * Xk[i];
		}

		Xk_ = Xk;
		Xk = Xk1;
		++k;
		if (k == 1)
			kEst = k_est(normC, normInf(Xk1 - Xk_), eps);
	} while (nevazka(A, B, Xk1) > eps);
	cout << "k = " << k << " <= " << kEst << endl;
	cout << "|| Ошибка || = " << normInf(Xk1 - Test2) << endl;

	return Xk1;
}

/****
Методы Зейделя и релаксации для трехдиагональных матриц
****/

// Обратная к нижнетругольной двухдиагонльной матрице (главная диагональ и
// диагональ ниже, заполняет строки не до конца, результирующая матрица будет 
// нижнетреугольной, не двухдиагональной).
// Использутеся для вычисления матрицы C в методах Зейделя и релаксации для
// трехдиагональных матриц
template<typename T>
matrix<T> inverseL2d(const vector<T>& a, const vector<T>& b) {
	size_t n = b.size();
	matrix<T> res;
	res.reserve(n);

	vector<T> row1{ 1 / b[0] };
	res.push_back(row1);
	for (size_t i = 1; i < n; ++i) {
		vector<T> row;
		row.reserve(i + 1);

		T resii = 1 / b[i];
		T coef = -a[i - 1] / b[i];

		for (size_t j = 0; j < i; ++j)
			row.push_back(coef * res[i - 1][j]);

		row.push_back(resii);
		res.push_back(row);
	}
	return res;
}

// Произведение нижнетреугольной и верхнетреугольной двухдиагональной матрицы матриц.
// Метод предполагает, что элементы под диагональю верхнетреугольной матрицы не
// заполнены (вместо них идут сразу диагольный элемент и т.д.)
// Использутеся для вычисления матрицы C в методах Зейделя и релаксации для
// трехдиагональных матриц.
// Ub - главаная диагональ
// Uc - диагональ над главной диагональю
template<typename T>
matrix<T> multLU2d(const matrix<T>& L, const vector<T>& Ub, const vector<T>& Uc) {
	size_t n = L.size();
	matrix<T> res;
	res.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row;
		row.reserve(n);

		// заполняем первый столбец
		T el = L[i][0] * Ub[0];
		row.push_back(el);

		// далее до диагонали включительно
		for (size_t j = 1; j <= i; ++j) {
			T el = L[i][j] * Ub[j] + L[i][j - 1] * Uc[j - 1];
			row.push_back(el);
		}

		// заполняем диагональ над диагональю
		if (i != n - 1) {
			T el = L[i][i] * Uc[i];
			row.push_back(el);
		}
		// дополняем нулями
		// делаем проверку, чтобы не было переполнения от size_t
		if (i < n - 2)
			for (size_t j = 0; j < n - (i + 2); ++j)
				row.push_back(0);
		res.push_back(row);
	}
	return res;
}

// Возвращает матрицу С для методов Зейделя (omega = 1) и релаксации для
// трехдиагональной матрицы
template<typename T>
matrix<T> relaxC3d(
	const vector<T>& a,
	const vector<T>& b,
	const vector<T>& c,
	T omega
) {
	vector<T> wL = a * omega;
	matrix<T> DwLinv = inverseL2d(wL, b);
	vector<T> D1w = b * (1 - omega);
	vector<T> wU = c * (-omega);
	return multLU2d(DwLinv, D1w, wU);
}

// Метод Зейделя для трехдиагональной матрицы
template<typename T>
vector<T> SeidelMethod3d(
	//const matrix<T>& A,
	const vector<T>& a,
	const vector<T>& b,
	const vector<T>& c,
	const vector<T>& B,
	const vector<T>& X0,
	T eps
) {
	int k = 0;
	size_t n = b.size();
	matrix<T> C = relaxC3d(a, b, c, T(1.0));
	//print_vec(C);
	T normC = normInfMatr(C);
	cout << endl << "||C|| = " << normC << endl;

	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1(n, 0);
	// норма невязки

	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;
	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}
	int kEst = 0;
	do
	{
		Xk1[0] = (B[0] - c[0] * Xk[1]) / b[0];
		for (size_t i = 1; i < n - 1; ++i) {
			Xk1[i] = (B[i] - a[i - 1] * Xk1[i - 1] - c[i] * Xk[i + 1]) / b[i];
		}
		Xk1[n - 1] = (B[n - 1] - a[n - 2] * Xk1[n - 2]) / b[n - 1];
		Xk_ = Xk;
		++k;
		Xk = Xk1;
		if (k == 1)
			kEst = k_est(normC, normInf(Xk1 - Xk_), eps);

	} while (stop(Xk_, Xk1, coef, eps));
	cout << "k = " << k << " <= " << kEst << endl;


	return Xk1;
}


// Метод релаксации для трехдиагональной матрицы
template<typename T>
vector<T> RelaxMethod3d(
	//const matrix<T>& A,
	const vector<T>& a,
	const vector<T>& b,
	const vector<T>& c,
	const vector<T>& B,
	const vector<T>& X0,
	T omega,
	T eps
) {
	int k = 0;
	size_t n = b.size();
	matrix<T> C = relaxC3d(a, b, c, omega);
	T normC = normInfMatr(C);
	cout << endl << "||C|| = " << normC << endl;


	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1(n, 0);
	// норма невязки

	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;
	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}
	int kEst;
	do
	{
		Xk1[0] = omega * (B[0] - c[0] * Xk[1]) / b[0] + (1 - omega) * Xk[0];
		for (size_t i = 1; i < n - 1; ++i) {
			Xk1[i] = omega * (B[i] - a[i - 1] * Xk1[i - 1] - c[i] * Xk[i + 1])
				/ b[i] + (1 - omega) * Xk[i];
		}
		Xk1[n - 1] = omega * (B[n - 1] - a[n - 2] * Xk1[n - 2]) / b[n - 1] + (1 - omega) * Xk[n - 1];

		Xk_ = Xk;
		++k;
		Xk = Xk1;
		if (k == 1)
			kEst = k_est(normC, normInf(Xk1 - Xk_), eps);
		// TODO критерий выхода с матрицей C
	} while (stop(Xk_, Xk1, coef, eps));
	cout << "k = " << k << " <= " << kEst << endl;

	return Xk1;
}

// Делает оценку сверху для модуля собственного значения матрицы
// Может быть полезным для настройки tau для метода простых итераций
template<typename T>
T GershgorinMax(const matrix<T>& A) {
	size_t n = A.size();
	T lambdaMax = 0;
	for (size_t i = 0; i < n; ++i) {
		T center = A[i][i];
		T radius = 0;
		for (size_t j = 0; j < n; ++j)
			radius += abs(A[i][j]);
		radius -= abs(center);
		T newLambda = abs(center) + radius;
		lambdaMax = max(lambdaMax, newLambda);
	}
	return lambdaMax;
}



// Раскладывает матрицу A = L + D + U, где L - нижнетруегольная, D - диагональная,
// U - верхнетреугольная
template<typename T>
tuple<matrix<T>, matrix<T>, matrix<T>> LDUdecomposition(const matrix<T>& A) {
	size_t n = A.size();

	matrix<T> L;
	L.reserve(n);
	L.push_back(vector<T>(n, 0));
	for (size_t i = 1; i < n - 1; ++i) {
		vector<T> row;
		row.reserve(n);
		for (size_t j = 0; j < i; ++j)
			row.push_back(A[i][j]);

		for (size_t j = 0; j < n - i; ++j)
			row.push_back(0);

		L.push_back(row);
	}
	L.push_back(A[n - 1]);
	L[n - 1][n - 1] = 0;

	matrix<T> D = zero<T>(n);
	for (size_t i = 0; i < n; ++i)
		D[i][i] = A[i][i];

	matrix<T> U;
	U.reserve(n);
	U.push_back(A[0]);
	U[0][0] = 0;
	for (size_t i = 1; i < n - 1; ++i) {
		vector<T> row;
		row.reserve(n);
		for (size_t j = 0; j <= i; ++j)
			row.push_back(0);

		for (size_t j = 0; j < n - i - 1; ++j)
			row.push_back(A[i][i + 1 + j]);

		U.push_back(row);
	}
	U.push_back(vector<T>(n, 0));

	return make_tuple(L, D, U);
}

// Выполнение алгоритмов для трехдиагональной матрицы 3-го теста
template<typename T>
void test3print(size_t N, T omega, T eps) {
	auto [d1, d2, d3, B] = test3d<T>(N);

	//print_vec(d1);
	//cout << endl;

	//print_vec(d2);
	//cout << endl;

	//print_vec(d3);
	//cout << endl;

	//print_vec(B);
	//cout << endl;

	auto [vec_A, vec_B] = test3<T>(N);
	vector<T> X0(200 + N, 5);

	// TODO для третьего теста последние несколько элементов отличаются от правильных
	cout << "Метод Зейделя для трехдиагональной матрицы:";
	//auto XSeidel = SeidelMethod3d(vec_A, vec_B, X0, eps);
	auto XSeidel = SeidelMethod3d(d1, d2, d3, vec_B, X0, eps);
	//print_vec(XSeidel);
	cout << endl;

	cout << "----------------------------------" << endl;
	cout << "Метод Релаксации для трехдиагональной матрицы:";
	//auto XRelax = RelaxMethod3d(vec_A, vec_B, X0, omega, eps);
	auto XRelax = RelaxMethod3d(d1, d2, d3, vec_B, X0, omega, eps);
	//print_vec(XRelax);
	cout << endl;
}



// Метод Зейделя для трехдиагональной матрицы
template<typename T>
vector<T> SeidelMethod3dopt(
	//const matrix<T>& A,
	const vector<T>& a,
	const vector<T>& b,
	const vector<T>& c,
	const vector<T>& B,
	const vector<T>& X0,
	T eps
) {
	int k = 0;
	size_t n = b.size();
	matrix<T> C = relaxC3d(a, b, c, T(1.0));
	//print_vec(C);
	T normC = normInfMatr(C);
	cout << endl << "||C|| = " << normC << endl;

	vector<T> Xk(X0);
	vector<T> Xk_;
	vector<T> Xk1(n, 0);
	// норма невязки

	// коэффициент для критерия выхода
	T coef = (1 - normC) / normC * eps;
	if (normC >= 1) {
		cout << endl << "Большая вероятность, что метод не сойдется. ||C|| = " << normC << endl;
		// из-за отрицательного коэффициента используем стандартный критерий, чтобы
		// хотя бы попытаться
		coef = eps;
	}
	int kEst = 0;
	do
	{
		if (k % 2 == 0) {
			Xk1[0] = (B[0] - c[0] * Xk[1]) / b[0];
			for (size_t i = 1; i < n - 1; ++i) {
				Xk1[i] = (B[i] - a[i - 1] * Xk1[i - 1] - c[i] * Xk[i + 1]) / b[i];
			}
			Xk1[n - 1] = (B[n - 1] - a[n - 2] * Xk1[n - 2]) / b[n - 1];
		}
		else
		{
			Xk1 = vector<T>(n, 0);
			Xk1[n - 1] = (B[n - 1] - a[n - 2] * Xk[n - 2]) / b[n - 1];
			for (size_t i = n - 2; i >= 1; i--) {
				Xk1[i] = (B[i] - a[i - 1] * Xk[i - 1] - c[i] * Xk1[i + 1]) / b[i];
			}
			Xk1[0] = (B[0] - c[0] * Xk1[1]) / b[0];
		}
		Xk_ = Xk;
		++k;
		Xk = Xk1;
		if (k == 1)
			kEst = k_est(normC, normInf(Xk1 - Xk_), eps);

	} while (normInf(Xk_ - Xk1) > eps);
	cout << "k = " << k << " <= " << kEst << endl;


	return Xk1;
}




template<typename T>
void example(size_t N, T omega, T eps) {
	vector<T> a1(N - 1, 6);
	vector<T> b1(N, 8);
	vector<T> c1(N - 1, 1);

	vector<T> a2(N - 1, 1);
	vector<T> b2(N, 8);
	vector<T> c2(N - 1, 6);
	vector<T> d;
	vector<T> X0(N, 1);

	for (size_t i = 0; i < N; i++) {
		d.push_back(i + 1);
	}
	auto XSeidel1 = SeidelMethod3d(a1, b1, c1, d, X0, eps);
	//print_vec(XSeidel1);

	auto XSeidel2 = SeidelMethod3d(a2, b2, c2, d, X0, eps);
	//print_vec(XSeidel2);

	/*auto XSeidel3 = RelaxMethod3d(a1, b1, c1, d, X0, omega, eps);
	print_vec(XSeidel3);

	auto XSeidel4 = RelaxMethod3d(a2, b2, c2, d, X0, omega, eps);
	print_vec(XSeidel4);*/

	auto XSeidel5 = SeidelMethod3dopt(a1, b1, c1, d, X0, eps);
	//print_vec(XSeidel5);

	auto XSeidel6 = SeidelMethod3dopt(a2, b2, c2, d, X0, eps);
	//print_vec(XSeidel6);
}




int main()
{
	setlocale(LC_ALL, "Russian");

	typedef double ddouble_t;

	//задаем вектора
	vector<vector<ddouble_t>> vec_A;// меняем в процессе работы
	vector<ddouble_t> vec_B;//меняем в процессе работы

	//считываем из файла размерность и матрицу
	//size_t N = out_of_file(vec_A, vec_B, "Вариант 1.txt");
	//size_t N = out_of_file(vec_A, vec_B, "Тест_2.txt");
	size_t N = out_of_file(vec_A, vec_B, "Тест 3.txt");
	//size_t N = out_of_file(vec_A, vec_B, "P_DAT1.txt");
	//size_t N = out_of_file(vec_A, vec_B, "DATA1.txt");

	cout << "Расширенная матрица:";
	print_AB(vec_A, vec_B);

	vector<ddouble_t> X0{ 0.0, 0.0 };
	//vector<ddouble_t> X0{ 5.0, 5.0, 5.0, 5.0 };
	//vector<ddouble_t> X0{ 500.0, 500.0, 500.0, 500.0 };
	ddouble_t eps = 1e-3;


	//cout << endl << "Метод простой итерации:";
	// Оценка для tau
	/*ddouble_t m = GershgorinMax(vec_A);
	cout << endl << "tau <= " << 2 / m;
	auto Xsimple = simpleIterationMethod(vec_A, vec_B, X0, 2 / m, eps);*/

	ddouble_t tau;

	tau = 2.0 / 183.0;
	cout << endl << "tau = " << tau << endl;
	auto Xsimple0 = simpleIterationMethod(vec_A, vec_B, X0, tau, eps);
	print_vec(Xsimple0);
	cout << "Невязка: " << nevazka(vec_A, vec_B, Xsimple0) << endl;
	cout << endl;

	/*tau = 0.0071518;
	cout << endl << "tau = " << tau << endl;
	auto Xsimple1 = simpleIterationMethod(vec_A, vec_B, X0, tau, eps);
	print_vec(Xsimple1);
	cout << "Невязка: " << nevazka(vec_A, vec_B, Xsimple1) << endl;
	cout << endl;


	tau = 0.001;
	cout << endl << "tau = " << tau << endl;
	auto Xsimple2 = simpleIterationMethod(vec_A, vec_B, X0, tau, eps);
	print_vec(Xsimple2);
	cout << "Невязка: " << nevazka(vec_A, vec_B, Xsimple2) << endl;
	cout << endl;


	tau = 0.0092558;
	cout << endl << "tau = " << tau << endl;
	auto Xsimple3 = simpleIterationMethod(vec_A, vec_B, X0, tau, eps);
	print_vec(Xsimple3);
	cout << "Невязка: " << nevazka(vec_A, vec_B, Xsimple3) << endl;
	cout << endl;*/

	/*tau = 0.02;
	cout << endl << "tau = " << tau << endl;
	auto Xsimple2 = simpleIterationMethod(vec_A, vec_B, X0, tau, eps);
	print_vec(Xsimple2);
	cout << "Невязка: " << nevazka(vec_A, vec_B, Xsimple2) << endl;
	cout << endl;*/


	//cout << "Метод Якоби:" << endl;
	//auto XJacobi = JacobiMethod(vec_A, vec_B, X0, eps);
	//print_vec(XJacobi);
	//cout << "Невязка: " << nevazka(vec_A, vec_B, XJacobi) << endl;
	//cout << endl;



	//cout << "Метод Зейделя:" << endl;
	//auto XSeidel = SeidelMethod(vec_A, vec_B, X0, eps);
	//print_vec(XSeidel);
	//cout << "Невязка: " << nevazka(vec_A, vec_B, XSeidel) << endl;
	//cout << endl;



	//ddouble_t omega1 = 1.1;
	//cout << "Метод релаксации, omega = " << omega1 << ":" << endl;
	//auto XRelax1 = RelaxMethod(vec_A, vec_B, X0, omega1, eps);
	//print_vec(XRelax1);
	//cout << "Невязка: " << nevazka(vec_A, vec_B, XRelax1) << endl;
	//cout << endl;


	//ddouble_t omega2 = 0.5;
	//cout << "Метод релаксации, omega = " << omega2 << ":" << endl;
	//auto XRelax2 = RelaxMethod(vec_A, vec_B, X0, omega2, eps);
	//print_vec(XRelax2);
	//cout << "Невязка: " << nevazka(vec_A, vec_B, XRelax2) << endl;
	//cout << endl;



	//test3print(15, 1.0, eps);
	//test3print(2, 0.5, eps);


	/*auto [L, D, U] = LDUdecomposition(vec_A);
	print_vec(vec_A);
	cout << endl;

	print_vec(L);
	cout << endl;

	print_vec(D);
	cout << endl;

	print_vec(U);
	cout << endl;*/

	cout << endl;

	//example(1000, 1.1, 0.0001);

}


