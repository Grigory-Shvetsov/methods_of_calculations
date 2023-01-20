#include<iostream>
#include<vector>
#include<iomanip>
#include<fstream>

using namespace std;

template<typename T>
using matrix = vector<vector<T>>;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

//считает сумму, чего надо вычесть
template<typename T>
T sum(size_t i, size_t N, const matrix<T>& vec_A, const vector<T>& vec_X)
{
	T sum = 0;
	int k = 0;
	for (size_t j = i + 1; j < N; j++)
	{
		sum += vec_A[i][j] * vec_X[k];
		k++;
	}
	return sum;
}

//обратный ход
// сложность O(n^2 / 2)
template<typename T>
void reverse_steps(int N, const matrix<T>& vec_A, const vector<T>& vec_B, vector<T>& vec_X)
{
	T x = vec_B[N - 1] / vec_A[N - 1][N - 1];
	vec_X.insert(vec_X.begin(), x);

	for (int i = N - 2; i > -1; i--) {
		x = (vec_B[i] - sum(i, N, vec_A, vec_X)) / vec_A[i][i];
		vec_X.insert(vec_X.begin(), x);
	}
}

//перестановка строк (частичный выбор)
template<typename T>
void main_str(size_t i, vector<vector<T>>& vec_A, vector<T>& vec_B)
{
	size_t k = i;
	// Ищем максимум
	T max = abs(vec_A[i][i]);
	for (size_t j = i + 1; j < vec_B.size(); j++) {
		// TODO стоит проверять на ноль (но что использовать в качетсве нуля для чисел с плавающей запятой)
		if (abs(vec_A[j][i]) >= max) {
			max = abs(vec_A[j][i]);
			k = j;
		}
	}

	// тут поменял
	if (max <= numeric_limits<T>::epsilon() * 10)
		throw "Вырожденная матрица";

	// Меняем строки
	swap(vec_A[k], vec_A[i]);
	swap(vec_B[k], vec_B[i]);
}


//Метод Гаусса
template<typename T>
vector<T> gauss(vector<vector<T>>& A, vector<T>& B) {
	// Поверим, что матрица A квадратная и имеет такой же размер, что и у B
	size_t N = A.size();
	vector<T> X;
	// тут поменял
	T c_ik;
	for (size_t k = 0; k < N - 1; ++k) {
		main_str(k, A, B);//строку поднимаем с max-элементом
		for (size_t i = k + 1; i < N; ++i) {
			//for (size_t i = k; i < N; ++i) {
			c_ik = A[i][k] / A[k][k];

			B[i] -= c_ik * B[k];
			for (size_t j = k + 1; j < N; ++j)
				A[i][j] -= c_ik * A[k][j];
		}
	}

	// При достаточном приближении к собственному числу матрица становится
	// похожей на вырожденную. Поэтому проверку лучше закомментировать
	// 
	// вроде бы это еще проверяет ранг (то, что ранг матрицы A равен рангу расширенной матрицы)
	//if (abs(A[N - 1][N - 1]) <= numeric_limits<T>::epsilon() * 10)
	//	throw "Вырожденная матрица";

	reverse_steps(N, A, B, X); // обратный ход
	return X;
}


// умножение матрицы и столбца
template<typename T>
vector<T> mult(const matrix<T>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}

// Скалярное произведение двух векторов
template<typename T>
T mult(const vector<T>& A1, const vector<T>& A2) {
	size_t n = A1.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		res += A1[i] * A2[i];
	return res;
}


// умножение матриц
template<typename T>
vector<vector<T>> mult(const matrix<T>& A1, const matrix<T>& A2) {
	size_t l = A1.size();
	size_t m = A2.size();
	size_t n = A2[0].size();

	vector<vector<T>> res;
	res.reserve(l);
	for (size_t i = 0; i < l; ++i) {
		vector<T> row;
		row.reserve(n); // оптимизируем память
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


// Создает единичную матрицу размером NxN
template<typename T>
vector<vector<T>> eig(size_t n) {
	vector<vector<T>> e;
	e.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		vector<T> row(n, 0);
		row[i] = 1;
		e.push_back(row);
	}
	return e;
}


// транспонирование квадратной матрицы (изменяет входную матрицу на транспорированую)
template<typename T>
void transpose(vector<vector<T>>& A) {
	size_t n = A.size();
	// тут была ошибка
	for (size_t i = 0; i < n - 1; ++i)
		for (size_t j = i + 1; j < n; ++j)
			swap(A[i][j], A[j][i]);
}



//Вывод матрицы А
template<typename T>
void print_A(const matrix<T>& vec_A)
{
	for (auto& i : vec_A) {
		for (auto j : i)
			cout << fixed << setprecision(4) << j << "\t";
		cout << endl;
	}
}


//Вывод вектора
template<typename T>
void print_vec(const vector<T>& vec)
{

	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << fixed << setprecision(8) << vec[i] << "\t";
	cout << endl << endl;
}

//считывание из файла
template<typename T>
size_t out_of_file(vector<vector<T>>& vec_A, string str)
{
	int N = 0;
	ifstream fout(str);
	fout >> N;
	T x_A;
	for (size_t i = 0; i < N; ++i) {
		vector<T> v_A;
		for (size_t j = 0; j < N; ++j) {
			fout >> x_A;
			v_A.push_back(x_A);
		}
		vec_A.push_back(v_A);
	}
	fout.close();
	return N;
}

// Умножение на матрицу поворота T_ij слева
// initIndex указывает на то, с какого столбца вычислять значения (нужно для матрицы R в QR-разложении)
template<typename T>
void multT(vector<vector<T>>& A, size_t n, size_t ti, size_t tj, T c, T s, size_t initIndex = 0) {
	//size_t n = A.size();
	// Копируем, потому что в цикле эти строки будут меняться
	vector<T> Ati(A[ti]);
	vector<T> Atj(A[tj]);

	for (size_t j = initIndex; j < n; ++j) {
		A[ti][j] = c * Ati[j] + s * Atj[j];
		A[tj][j] = -s * Ati[j] + c * Atj[j];
	}
}


// раскладывает матрицу A в QR-разложение, где Q - ортогональная матрица, R - верхнедиагональная (записывается в A)
template<typename T>
void qr(vector<vector<T>>& A, vector<vector<T>>& Q)
{
	// для этой функции матрица A должна быть квадратной
	size_t n = A.size();
	Q = eig<T>(n);

	// матрица состоит из чисел c и s
	for (size_t i = 0; i < n - 1; ++i)
	{

		for (size_t j = i + 1; j < n; ++j) {
			T Aii = A[i][i];
			T Aji = A[j][i];

			// знаментель
			T denom = sqrt(Aii * Aii + Aji * Aji);

			if (denom <= numeric_limits<T>::epsilon() * 10)
				continue;

			T c_ij = Aii / denom;
			T s_ij = Aji / denom;

			// тут поменял
			multT(Q, n, i, j, c_ij, s_ij);
			multT(A, n, i, j, c_ij, s_ij, i);
		}
	}
	// тут поменял
	// проверка на вырожденность
	for (size_t i = 0; i < n; ++i)
		if (abs(A[i][i]) <= numeric_limits<T>::epsilon() * 10)
			throw "Вырожденная матрица";
}

// норма (евклидова) разности двух векторов
template<typename T>
T norm(const vector<T>& x1, const vector<T>& x2) {
	T res = 0;
	size_t n = x1.size();
	for (size_t i = 0; i < n; ++i)
		res += (x1[i] - x2[i]) * (x1[i] - x2[i]);
	return sqrt(res);
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

// нормирует вектор
template<typename T>
void normalize(vector<T>& x) {
	size_t n = x.size();
	T normX = norm(x);
	for (size_t i = 0; i < n; ++i)
		x[i] /= normX;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

template<typename T>
vector<T> eigenvector(const matrix<T> A, T eigVal, T eps) {
	size_t n = A.size();

	// начальное приближение
	vector<T> xk(n, 0);
	xk[0] = 1;

	//print_vec(xk);

	// матрица для решения СЛАУ
	matrix<T> M(A);
	for (size_t i = 0; i < n; ++i)
		M[i][i] -= eigVal;

	vector<T> xk1;

	T coef;
	do {
		matrix<T> MCopy(M);
		vector<T> xkCopy(xk);
		xk1 = gauss(MCopy, xkCopy);
		normalize(xk1);
		coef = norm(xk1, xk);

		// Может возникнуть такая ситуация, что начнется чередоваться один и
		// тот же вектор, просто домноженный на -1
		for (size_t i = 0; i < n; ++i)
			xk[i] *= -1;
		coef = min(coef, norm(xk1, xk));

		xk = xk1;

		//print_vec(xk);
	} while (coef >= eps);

	return xk;
}

// Получает матрицу Хессенберга
//если матрица симметрична, можно оптимизировать
template<typename T>
void Hessenberg(matrix<T>& A)
{
	// для этой функции матрица A должна быть квадратной
	size_t n = A.size();

	for (size_t i = 1; i < n - 1; ++i)
	{
		for (size_t j = i + 1; j < n; ++j) {
			T c = A[i][i - 1];
			T s = A[j][i - 1];

			// знаментель
			T denom = sqrt(c * c + s * s);

			//if (denom <= numeric_limits<T>::epsilon() * 10)
			//	continue;
			T alpha = c / denom;
			T beta = s / denom;

			//умножение матрицы А слева на T_kl
			//for (size_t k = 0; k < n; k++) {
			for (size_t k = i - 1; k < n; k++) {

				// Копируем, потому что эти элементы будут меняться
				T Aik = A[i][k];
				T Ajk = A[j][k];

				A[i][k] = alpha * Aik + beta * Ajk;
				A[j][k] = -beta * Aik + alpha * Ajk;
			}


			//умножение справа на (T_kl)^T получившейся матрицы
			for (size_t k = 0; k < n; k++) {
				//сохраняем значения, они меняются
				T Aki(A[k][i]);
				T Akj(A[k][j]);

				A[k][i] = alpha * Aki + beta * Akj;
				A[k][j] = -beta * Aki + alpha * Akj;
			}

		}
	}

}



// Ищет RQ для матрицы Хессенберга. Меняет входную матрицу
// n - какую размерность для матрицы A следует иметь в виду
// в процессе поиска собственных чисел n будет уменьшатся
template<typename T>
void qrHessenberg(matrix<T>& A, size_t n) {
	for (size_t i = 0; i < n - 1; ++i) {
		size_t j = i + 1;

		T c = A[i][i];
		T s = A[j][i];

		// знаментель
		T denom = sqrt(c * c + s * s);
		T alpha = c / denom;
		T beta = s / denom;

		// Копируем, потому что в цикле эти строки будут меняться
		vector<T> Aik(A[i]);
		vector<T> Ajk(A[j]);
		//умножение матрицы А слева на T_kl
		//for (size_t k = i; k < n; k++) {
		// Здесь нужно учесть, что 
		for (size_t k = min(i - 1, i); k < n; k++) {
			A[i][k] = alpha * Aik[k] + beta * Ajk[k];
			A[j][k] = -beta * Aik[k] + alpha * Ajk[k];
		}

		//умножение справа на (T_kl)^T получившейся матрицы
		// TODO мб можно сократить количество операций поменяв n на ...
		for (size_t k = 0; k < n; k++) {
			//сохраняем значения, они меняются
			T aki(A[k][i]);
			T akj(A[k][j]);

			A[k][i] = alpha * A[k][i] + beta * A[k][j];
			A[k][j] = -beta * aki + alpha * akj;
		}
	}
}



// Поиск собственных чисел методом QR-разложения
// Возвращает список собственных чисел
// TODO считать количество операций
template<typename T>
vector<T> eigValQR(const matrix<T>& A, T eps) {
	matrix<T> H(A);
	Hessenberg(H);
	size_t n = H.size();
	//cout << "Вывод матрицы Хессенберга перед QR-разложением" << endl;
	//print_A(H);
	//cout << endl;

	vector<T> res;
	res.reserve(n);

	size_t iterCount = 0;
	size_t operCount = 0;

	while (n != 1) {
		while (abs(H[n - 1][n - 2]) >= eps) {
			//выбираем шаг сдвига
			T sigma = H[n - 1][n - 1];
			//sigma = 0.0;
			//выполняем сдвиг (-сигма)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] -= sigma;
			}
			//раскладываем матрицу на A_k+1 = Q * A_k * (Q)^T 
			qrHessenberg(H, n);
			operCount += 6 * n * n + 7 * n - 13;
			//выполняем обратный сдвиг (+сигма)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] += sigma;
			}

			/*print_A(H);
			cout << endl*/;

			iterCount++;
		}
		res.push_back(H[n - 1][n - 1]);
		--n;
	}
	res.push_back(H[0][0]); //lambda 1 (без сдвига)

	cout << "Кол-во итераций: " << iterCount << endl;
	cout << "Кол-во операций: " << operCount << endl;
	return res;
}



// Ищет RQ для иходной матрицы (со/без сдвига). Меняет входную матрицу
// n - какую размерность для матрицы A следует иметь в виду
// в процессе поиска собственных чисел n будет уменьшатся
template<typename T>
void qrHessenbergSimple(matrix<T>& A, size_t n) {
	matrix<T> Q = eig<T>(n);
	for (size_t i = 0; i < n - 1; ++i) {
		for (size_t j = i + 1; j < n; ++j) {
			T Aii = A[i][i];
			T Aji = A[j][i];

			// знаментель
			T denom = sqrt(Aii * Aii + Aji * Aji);

			if (denom <= numeric_limits<T>::epsilon() * 10)
				continue;

			T c_ij = Aii / denom;
			T s_ij = Aji / denom;

			// тут поменял
			multT(Q, n, i, j, c_ij, s_ij);
			multT(A, n, i, j, c_ij, s_ij, i);
		}
	}
	transpose(Q);
	A = mult(A, Q);

}

template<typename T>
vector<T> eigValQRSimple(const matrix<T>& A, T eps) {
	matrix<T> H(A);
	size_t n = H.size();

	vector<T> res;
	res.reserve(n);

	size_t iterCount = 0;
	size_t operCount = 0;

	while (n != 1) {
		while (abs(H[n - 1][n - 2]) >= eps) {
			//выбираем шаг сдвига
			T sigma = H[n - 1][n - 1];
			//sigma = 0.0;
			//выполняем сдвиг (-сигма)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] -= sigma;
			}
			//раскладываем матрицу на A_k+1 = Q * A_k * (Q)^T 
			qrHessenbergSimple(H, n);

			operCount += 10.0 * n * n*n/3.0 + n * n/2.0 - 23*n/6.0;
			//выполняем обратный сдвиг (+сигма)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] += sigma;
			}

			/*print_A(H);
			cout << endl;*/

			iterCount++;
		}
		res.push_back(H[n - 1][n - 1]);
		--n;
	}
	res.push_back(H[0][0]); //lambda 1 (без сдвига)

	cout << "Кол-во итераций: " << iterCount << endl;
	cout << "Кол-во операций: " << operCount << endl;
	return res;
}


template<typename T>
vector<pair<T, vector<T>>> eigValVecQRSimple(const matrix<T>& A, T eps) {
	vector<T> eigVals = eigValQRSimple(A, eps);
	//reverse(eigVals.begin(), eigVals.end());

	vector<pair<T, vector<T>>> res;
	res.reserve(eigVals.size());
	for (size_t i = 0; i < eigVals.size(); ++i) {
		//cout << fixed << setprecision(8) << "lambda = " << eigVals[i] << ":" << endl;

		auto eigVec = eigenvector(A, eigVals[i], eps);
		res.push_back(make_pair(eigVals[i], eigVec));
	}

	return res;
}




// Поиск собственных чисел и собственных векторов с использованием QR-разложения
// и метода обратных операций.
// Возвращает список пар "Собственное число" - "Соответствующий ему собственный
// вектор".
template<typename T>
vector<pair<T, vector<T>>> eigValVecQR(const matrix<T>& A, T eps) {
	vector<T> eigVals = eigValQR(A, eps);
	//reverse(eigVals.begin(), eigVals.end());

	vector<pair<T, vector<T>>> res;
	res.reserve(eigVals.size());
	for (size_t i = 0; i < eigVals.size(); ++i) {
		//cout << fixed << setprecision(8) << "lambda = " << eigVals[i] << ":" << endl;

		auto eigVec = eigenvector(A, eigVals[i], eps);
		res.push_back(make_pair(eigVals[i], eigVec));
	}

	return res;
}

// генерирует начальное приближение к собственному вектору случайными числами
template<typename T>
vector<T> generateX0(size_t n) {
	vector<T> res;
	res.reserve(n);
	for (size_t i = 0; i < n; ++i)
		res.push_back(rand());
	normalize(res);
	return res;
}

// Поиск собственных чисел и векторов при помощи соотношения Рэлея
// Возвращает список пар "Собств. число"-"Собств. вектор"
template<typename T>
vector<pair<T, vector<T>>> RayleighRelation(const matrix<T>& A, T eps) {
	size_t n = A.size();
	vector<pair<T, vector<T>>> res;
	res.reserve(n);

	matrix<T> A1(A);

	for (size_t i = 0; i < n; ++i) {
		// начальное приближение собственного вектора
		vector<T> xk = generateX0<T>(n);
		// для первого соб. вектора пропускается
		// для следующих, чтобы не получить одного и того же соб. вектора
		for (size_t j = 0; j < res.size(); ++j) {
			vector<T> ej = res[j].second;
			T f = mult(xk, ej);

			for (size_t k = 0; k < n; ++k)
				xk[k] -= f * ej[k];
		}
		normalize(xk);
		T lamk;
		vector<T> xk1;

		T coef;
		do {
			lamk = mult(mult(A1, xk), xk);

			// матрица для решения СЛАУ
			matrix<T> M(A1);
			for (size_t i = 0; i < n; ++i)
				M[i][i] -= lamk;

			vector<T> xkCopy(xk);
			xk1 = gauss(M, xkCopy);

			// удаляем из xk1 компоненты уже известных собственных векторов
			// с помощью коэффициентов ряда фурье
			for (size_t j = 0; j < res.size(); ++j) {
				vector<T> ej = res[j].second;
				T f = mult(xk1, ej);
				for (size_t k = 0; k < n; ++k)
					xk1[k] -= f * ej[k];
			}
			////////////////

			normalize(xk1);

			//критерий остановы норма разности последних приближений с учетом одного и
			//того же вектора с разными знаками
			coef = norm(xk1, xk);

			for (size_t i = 0; i < n; ++i)
				xk[i] *= -1;
			coef = min(coef, norm(xk1, xk));

			xk = xk1;
		} while (coef >= eps);

		lamk = mult(mult(A1, xk), xk);
		res.push_back(make_pair(lamk, xk));
	}
	return res;
}

// Проверяет, что число eigVal и вектор eigVec действительно являются
// собственными. Проверяет, что (A - eigVal * E) * eigVec = 0.
template<typename T>
bool checkEigValAndVec(const matrix<T>& A, T eigVal, const vector<T>& eigVec, T eps) {
	size_t n = A.size();
	matrix<T> AlamE(A);
	for (size_t i = 0; i < n; ++i) {
		AlamE[i][i] -= eigVal;
	}

	vector<T> mul = mult(AlamE, eigVec);
	for (size_t i = 0; i < mul.size(); ++i)
		if (abs(mul[i]) >= eps)
			return false;
	return true;
}



/////////////////////////////////////////////////////////////////////
// "Неправильное" приведение к матрице Хессенберга
// Пытается привести матрицу к верхнетреугольному виду
// (но ничего не получится)
// Для 4 контрольного вопроса
template<typename T>
void HessenbergIncorr(const matrix<T>& M)
{
	matrix<T> A(M);

	// для этой функции матрица A должна быть квадратной
	size_t n = A.size();

	//for (size_t i = 1; i < n - 1; ++i)
	for (size_t i = 0; i < n - 1; ++i)
	{
		for (size_t j = i + 1; j < n; ++j) {
			T c = A[i][i];
			T s = A[j][i];

			// знаментель
			T denom = sqrt(c * c + s * s);

			T alpha = c / denom;
			T beta = s / denom;

			//умножение матрицы А слева на T_kl
			for (size_t k = 0; k < n; k++) {
				// нули пробадают, поэтому ничего нельзя оптимизировать
				//for (size_t k = i; k < n; k++) {

					// Копируем, потому что эти элементы будут меняться
				T Aik = A[i][k];
				T Ajk = A[j][k];

				A[i][k] = alpha * Aik + beta * Ajk;
				A[j][k] = -beta * Aik + alpha * Ajk;
			}


			//умножение справа на (T_kl)^T получившейся матрицы
			for (size_t k = 0; k < n; k++) {
				//сохраняем значения, они меняются
				T Aki(A[k][i]);
				T Akj(A[k][j]);

				A[k][i] = alpha * Aki + beta * Akj;
				A[k][j] = -beta * Aki + alpha * Akj;
			}
			//print_A(A);

		}
	}
	cout << "\"Неправильный\" Хессенберг:" << endl;
	print_A(A);
	cout << "\"Правильный\" Хессенберг:" << endl;
	matrix<T> A2(M);
	Hessenberg(A2);
	print_A(A2);
}



// метод обратных итераций с выбором начального приближения
template<typename T>
vector<T> eigenvector(const matrix<T> A, T eigVal, const vector<T>& x0, T eps) {
	size_t n = A.size();

	// начальное приближение
	vector<T> xk = x0;

	// матрица для решения СЛАУ
	matrix<T> M(A);
	for (size_t i = 0; i < n; ++i)
		M[i][i] -= eigVal;

	vector<T> xk1;

	T coef;
	do {
		matrix<T> MCopy(M);
		vector<T> xkCopy(xk);
		xk1 = gauss(MCopy, xkCopy);
		normalize(xk1);
		coef = norm(xk1, xk);

		// Может возникнуть такая ситуация, что начнется чередоваться один и
		// тот же вектор, просто домноженный на -1
		for (size_t i = 0; i < n; ++i)
			xk[i] *= -1;
		coef = min(coef, norm(xk1, xk));

		xk = xk1;
	} while (coef >= eps);

	return xk;
}


// Поиск одного собственного чисел и вектора при помощи соотношения Рэлея
// Возвращает пару "Собств. число"-"Собств. вектор"
// Принимает на вход начальное приближение x0
template<typename T>
pair<T, vector<T>> RayleighRelationSimple(const matrix<T>& A, const vector<T>& x0, T eps) {
	size_t n = A.size();

	matrix<T> A1(A);


	// начальное приближение
	vector<T> xk = x0;
	normalize(xk);

	T lamk;

	vector<T> xk1;

	T coef;
	do {
		lamk = mult(mult(A1, xk), xk);

		// матрица для решения СЛАУ
		matrix<T> M(A1);
		for (size_t i = 0; i < n; ++i)
			M[i][i] -= lamk;

		vector<T> xkCopy(xk);
		xk1 = gauss(M, xkCopy);

		normalize(xk1);
		coef = norm(xk1, xk);

		for (size_t i = 0; i < n; ++i)
			xk[i] *= -1;
		coef = min(coef, norm(xk1, xk));

		xk = xk1;
		//print_vec(xk);
	} while (coef >= eps);

	lamk = mult(mult(A1, xk), xk);
	return make_pair(lamk, xk);
}

// Пример к 6 контрольному вопросу
void question6() {
	matrix<double> A;
	size_t N = out_of_file(A, "EIGEN1.TXT");

	double eps = 1e-3;

	double lam1 = 57.165516324888;
	// собственный вектор с собственным числом 57.165516324888
	vector<double> x1{ -0.063563494887, 0.443928222700, 0.888001562853, -0.101688935383 };

	double lam2 = 22.359615081837;
	// собственный вектор с собственным числом 22.359615081837
	vector<double> x2{ 0.777968847720, 0.582505207755, -0.229455791959, 0.052935757783 };

	cout << "Метод обратной итерации:";
	// сойдется к x1 (из-за lam1)
	auto eInv = eigenvector(A, lam1, x2, eps);
	print_vec(eInv);
	cout << "Соотношение Рэлея:";
	// сойдется к x2 (из-за начального приближения)
	auto lamERayleigh = RayleighRelationSimple(A, x2, eps);
	auto eRayleigh = lamERayleigh.second;
	print_vec(eRayleigh);
}



//метод QR-разложения
//R = P^(-1) * A * P
//P - невырожденная
//
//
//QR - разложение матрицы требует (O(n3)) операций
//значит нужно ускорить QR-разложение -> приводим матрицу А к форме Хессенберга

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << fixed;
	cout.precision(10);


	typedef double MyType;
	MyType eps = 1e-3;
	matrix<MyType> A;
	matrix<MyType> Q;

	//читаем из файла
	//size_t N = out_of_file(A, "Тест 1.txt");
	//size_t N = out_of_file(A, "Тест 2.txt"); //несимметричная матрица
	//size_t N = out_of_file(A, "EIGEN1.TXT");
	size_t N = out_of_file(A, "EIGEN20.TXT");



	//matrix<MyType> A_H(A);
	//cout << "Матрица А:" << endl;
	//print_A(A);
	//cout << endl << "Поиск матрицы Хессенберга" << endl;
	//Hessenberg(A_H);
	//cout << endl;
	//qrHessenberg(A_H, A_H.size());
	//print_A(A_H);

	//transpose(A);
	////QR-разложение
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\tQR - разложение" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	auto eigsQR = eigValVecQR(A, eps);
	for (auto& eig : eigsQR) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "Собственное число " << fixed << setprecision(10) << eigVal << endl << "Собственный вектор:";
		print_vec(eigVec);
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}


	/*QR-разложение (без матрицы Хессенберга)*/
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\tQR - разложение (без матрицы Хессенберга)" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	auto eigsQRSimple = eigValVecQRSimple(A, eps);
	for (auto& eig : eigsQRSimple) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "Собственное число " << fixed << setprecision(10) << eigVal << endl << "Собственный вектор:";
		print_vec(eigVec);
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}




	//соотношение Рэлея
	cout << endl << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\tСоотношение Рэлея" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	//srand(time(NULL));
	auto eigs = RayleighRelation(A, eps);
	for (auto& eig : eigs) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "Собственное число " << fixed << setprecision(10) << eigVal << endl << "Собственный вектор:";
		print_vec(eigVec);
		cout << endl;
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}


	question6();
}



//Ответ на 4 вопрос, как я понял, неправильный.Ответ вроде должен быть,
//что метод вращений мы не приведем к треугольному виду.По крайней мере, 
//когда я высказал эту идею, Гусев сказал, что идея хорошая.
//Для 6 вопроса сделать пример.
//Ещё в таблицу в результатах добавить количество операций и добавить для
//QR - алгоритма без приведения форме Хессенберга. И нужно показать, как 
//изменятся собственные числа и вектора для транспонированной матрицы
//(для несимметричной матрицы, конечно).Ну и привести пример
//Ещё когда Гусев спрашивал про собственные вектора транспонированной матрицы, 
//он зачем - то начал спрашивать про сопряженное пространство и теорему Рисса 
//о представлении функционала.Но зачем, я не понял