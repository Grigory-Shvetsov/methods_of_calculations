#include<iostream>
#include<iomanip>
#include<vector>
#include<string>
#include<fstream>
#include<algorithm>

using namespace std;

// косвенная индексация!!!

//Вывод матрицы А
template<typename T>
void print_A(const vector<vector<T>>& vec_A)
{
	cout << endl;
	for (auto& i : vec_A) {
		for (auto j : i)
			cout << j << "\t";
		cout << endl;
	}
}

template<typename T>
using matrix = vector<vector<T>>;

//Вывод расширенной матрицы
template<typename T>
void print_AB(const vector<vector<T>>& vec_A, const vector<T>& vec_B)
{
	cout << endl;
	for (size_t i = 0; i < vec_B.size(); i++) {
		for (size_t j = 0; j < vec_B.size(); j++)
			cout << fixed << setprecision(7) << vec_A[i][j] << "\t";
		cout << fixed << setprecision(7) << vec_B[i] << endl;
	}
}



//Вывод вектора
template<typename T>
void print_vec(const vector<T>& vec)
{
	cout << endl;
	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << fixed << setprecision(8) << vec[i] << "\t";
	cout << endl;
}


//считывание из файла
template<typename T>
size_t out_of_file(vector<vector<T>>& vec_A, vector<T>& vec_B, string str)
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


//считает сумму, чего надо вычесть
template<typename T>
T sum(size_t i, size_t N, const vector<vector<T>>& vec_A, const vector<T>& vec_X)
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
template<typename T>
void reverse_steps(int N, const vector<vector<T>>& vec_A, const vector<T>& vec_B, vector<T>& vec_X)
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

	// вроде бы это еще проверяет ранг (то, что ранг матрицы A равен рангу расширенной матрицы)
	if (abs(A[N - 1][N - 1]) <= numeric_limits<T>::epsilon() * 10)
		throw "Вырожденная матрица";

	reverse_steps(N, A, B, X); // обратный ход
	return X;
}



// умножение матрицы и столбца
template<typename T>
vector<T> mult(const vector<vector<T>>& A1, const vector<T>& A2) {
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



//Правая часть
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


//невязка
template<typename T>
vector<T> nevazka(const vector<vector<T>>& vec_A_source, const vector<T>& vec_B_source, const vector<T>& vec_X)
{
	vector<T> vec_AX = right_part(vec_A_source, vec_X);
	vector<T> vec;
	T bb1 = 0;
	for (size_t i = 0; i < vec_AX.size(); i++) {
		bb1 = vec_AX[i] - vec_B_source[i];
		vec.push_back(bb1);
	}
	return vec;
}


// умножение матриц
template<typename T>
vector<vector<T>> mult(vector<vector<T>>& A1, vector<vector<T>>& A2) {
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



// Умножение на матрицу поворота T_ij слева
// initIndex указывает на то, с какого столбца вычислять значения (нужно для матрицы R в QR-разложении)
template<typename T>
void multT(vector<vector<T>>& A, size_t ti, size_t tj, T c, T s, size_t initIndex = 0) {
	size_t n = A.size();
	// Копируем, потому что в цикле эти строки будут меняться
	vector<T> Ati(A[ti]);
	vector<T> Atj(A[tj]);

	for (size_t j = initIndex; j < n; ++j) {
		A[ti][j] = c * Ati[j] + s * Atj[j];
		A[tj][j] = -s * Ati[j] + c * Atj[j];
	}
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
			multT(Q, i, j, c_ij, s_ij);
			multT(A, i, j, c_ij, s_ij, i);
		}
	}
	// тут поменял
	// проверка на вырожденность
	for (size_t i = 0; i < n; ++i)
		if (abs(A[i][i]) <= numeric_limits<T>::epsilon() * 10)
			throw "Вырожденная матрица";
}



template<typename T>
// тут поменял
vector<T> fun_qr(vector<vector<T>>& A, vector<vector<T>>& Q, const vector<vector<T>>& vec_A_source, const vector<T>& vec_B_source)
{
	qr(A, Q);
	vector<T> vec_X2;
	reverse_steps(vec_B_source.size(), A, mult(Q, vec_B_source), vec_X2);
	/*transpose(Q);*/
	return vec_X2;
}




// норма p=2 квадратной матрицы
template <typename T>
T normMatr(const vector<vector<T>>& A)
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
T norm1Matr(const vector<vector<T>>& A)
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
T normInfMatr(const vector<vector<T>>& A)
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


//возмущение
template<typename T>
vector<T> outrage(vector<T>& vec, int k)
{
	vector<T> deltab(vec.size());
	switch (k) {
	case 2:
		for (size_t i = 0; i < vec.size(); i++) {
			vec[i] += 0.03;
			deltab[i] = 0.03;
		}
		break;
	case 1:
		for (size_t i = 0; i < vec.size(); i++) {
			vec[i] += 0.02;
			deltab[i] = 0.02;
		}
		break;
	case 0:
		for (size_t i = 0; i < vec.size(); i++) {
			vec[i] += 0.01;
			deltab[i] += 0.01;
		}
		break;
	default:
		throw "Некорректное значение k";
	}
	return deltab;
}




// находит обратную к квадратной матрице
template<typename T>
vector<vector<T>> inverse(const vector<vector<T>>& A)
{
	vector<vector<T>> AA(A);
	size_t n = AA.size();
	vector<vector<T>> Ainv;
	Ainv.reserve(n);
	vector<vector<T>> Q;
	qr(AA, Q);
	transpose(Q);
	for (size_t i = 0; i < n; ++i) {
		// столбец единичной матрицы
		//vector<T> e(n, 0);
		//e[i] = 1;
		// TODO умножение на столбец (0, 0, ..., 1, ..., 0) можно сделать проще
		//auto B = mult(Q, e);

		// тут поменял
		vector<T> B(n);
		for (size_t j = 0; j < n; ++j)
			B[j] = Q[j][i];

		vector<T> X;
		reverse_steps(n, AA, B, X);
		Ainv.push_back(X);
	}
	transpose(Ainv);
	return Ainv;
}



// обусловленность матрицы (меняет исходную матрицу)
template<typename T>
T cond1(vector<vector<T>>& A) {
	T norm1A = norm1Matr(A);

	auto Ainv = inverse(A);
	T norm1Ainv = norm1Matr(Ainv);

	return norm1A * norm1Ainv;
}


// обусловленность матрицы (меняет исходную матрицу)
template<typename T>
T condInf(vector<vector<T>>& A) {
	T norm1A = normInfMatr(A);

	auto Ainv = inverse(A);
	T norm1Ainv = normInfMatr(Ainv);

	return norm1A * norm1Ainv;
}


//оценка снизу (р = 1)
template<typename T>
T mark_down1(const vector<vector<T>>& A, const vector<T>& B, const vector<T>& X)
{
	vector<T> deltab, deltax;
	vector<vector<T>> vec_deltab(3), vec_deltax(3);

	//посчитали дельта-б и дельта-х
	for (size_t i = 0; i < 3; i++) {
		vector<T> B1(B);
		vec_deltab[i] = outrage(B1, i);
		//vec_deltax[i] = inverse(A) * vec_deltab[i];
		vec_deltax[i] = mult(inverse(A), vec_deltab[i]);
	}
	//посчитали отношения дельта-б и дельта-х
	for (size_t i = 0; i < 3; i++) {
		deltab.push_back(norm1(vec_deltab[i]) / norm1(B));
		deltax.push_back(norm1(vec_deltax[i]) / norm1(X));
	}

	vector<T> vec_deltax_deltab(3);
	//тут находим отношение дельта-x и дельта-b
	for (size_t i = 0; i < 3; i++) {
		vec_deltax_deltab[i] = deltax[i] / deltab[i];
	}
	//возвращаем максимум из трех возмущений(в книжке написано несколько раз и выбрать максимум)
	//return max(vec_deltax_deltab);
	return *max_element(vec_deltax_deltab.begin(), vec_deltax_deltab.end());
}



//оценка снизу (р = inf)
template<typename T>
T mark_downInf(const vector<vector<T>>& A, const vector<T>& B, const vector<T>& X)
{
	vector<T> deltab, deltax;
	vector<vector<T>> vec_deltab(3), vec_deltax(3);

	//посчитали дельта-б и дельта-х
	for (size_t i = 0; i < 3; i++) {
		vector<T> B1(B);
		vec_deltab[i] = outrage(B1, i);
		//vec_deltax[i] = inverse(A) * vec_deltab[i];
		vec_deltax[i] = mult(inverse(A), vec_deltab[i]);
	}
	//посчитали отношения дельта-б и дельта-х
	for (size_t i = 0; i < 3; i++) {
		deltab.push_back(normInf(vec_deltab[i]) / normInf(B));
		deltax.push_back(normInf(vec_deltax[i]) / normInf(X));
	}

	vector<T> vec_deltax_deltab(3);
	//тут находим отношение дельта-x и дельта-b
	for (size_t i = 0; i < 3; i++) {
		vec_deltax_deltab[i] = deltax[i] / deltab[i];
	}
	//возвращаем максимум из трех возмущений(в книжке написано несколько раз и выбрать максимум)
	//return max(vec_deltax_deltab);
	return *max_element(vec_deltax_deltab.begin(), vec_deltax_deltab.end());
}



// TODO сделать:
// 1. Оценка числа обусловленности

// обратная матрица:
// 1. находим QR-разложение, берем от них обратные
// 2. Решаем n СЛАУ
// 3. Соединяем их

int main()
{
	setlocale(LC_ALL, "Russian");

	typedef double ddouble_t;


	//задаем вектора
	vector<vector<ddouble_t>> vec_A;// меняем в процессе работы
	vector<ddouble_t> vec_B;//меняем в процессе работы

	//считываем из файла размерность и матрицу
	//size_t N = out_of_file(vec_A, vec_B, "Вариант 1.txt");
	size_t N = out_of_file(vec_A, vec_B, "Вариант 1.txt");

	// копируем матрицу и вектор, потому что метод их изменяет
	vector<ddouble_t> vec_B_source(vec_B);
	vector<vector<ddouble_t>> vec_A_source(vec_A);

	//нужны для функции fun_outrage - находит решение с возмущением
	vector<double_t> vec_B_outrage(vec_B);
	vector<vector<ddouble_t>> vec_A_outrage(vec_A);


	cout << "Расширенная матрица:";
	print_AB(vec_A, vec_B);


	//метод Гаусса (меняет матрицу и правую часть)
	auto vec_X1 = gauss(vec_A, vec_B);
	cout << endl << "Решение методом Гаусса:";
	print_vec(vec_X1); //решение


	//невязка p=1
	cout << endl << "Невязка р=1: ";
	cout << (norm1(nevazka(vec_A_source, vec_B_source, vec_X1)));

	//невязка p=2
	cout << endl << "Невязка р=2: ";
	cout << (norm(nevazka(vec_A_source, vec_B_source, vec_X1)));

	//невязка p=inf
	cout << endl << "Невязка р=inf: ";
	cout << (normInf(nevazka(vec_A_source, vec_B_source, vec_X1)));



	//вызов QR-разложения и вывод матриц Q и R
	vector<vector<ddouble_t>> Q;
	vector<vector<ddouble_t>> R = vec_A_source;
	auto vec_X2 = fun_qr(R, Q, vec_A_source, vec_B_source);
	cout << endl << "Решение QR-методом:";
	print_vec(vec_X2); //решение

	//матрицы Q и R
	cout << endl << "Матрица Q:";
	print_A(Q);
	cout << endl << "Матрица R:";
	print_A(R);

	//невязка p=1
	cout << endl << "Невязка р=1: ";
	cout << (norm1(nevazka(vec_A_source, vec_B_source, vec_X2)));

	//невязка p=2
	cout << endl << "Невязка р=2: ";
	cout << (norm(nevazka(vec_A_source, vec_B_source, vec_X2)));

	//невязка p=inf
	cout << endl << "Невязка р=inf: ";
	cout << (normInf(nevazka(vec_A_source, vec_B_source, vec_X2)));


	cout << endl;



	//обратная матрица
	/*vector<vector<ddouble_t>> A_inv(vec_A_source);
	vector<vector<ddouble_t>> vec_A_inv = inverse(A_inv);
	cout << "Обратная матрица:";
	print_A(vec_A_inv);*/

	//проверка обратной матрицы
	/*print_A(vec_A_source);
	print_A(mult(vec_A_inv,vec_A_source));*/


	ddouble_t cond_1 = cond1(vec_A_source);
	cout << "Обусловленность (p = 1): " << cond_1 << endl;

	//ddouble_t cond_1 = cond1(eig<double>(4));
	//cout << "Обусловленность (): " << cond_1 << endl;

	/*ddouble_t cond_Inf = condInf(vec_A_source);
	cout << "Обусловленность (p = inf): " << cond_Inf << endl;*/



	//оценка сверху
	//если знаем разложение матрицы A=A1*A2
	//оценка сверху cond A <= cond A1 * cond A2
	/*ddouble_t markup1 = cond1(Q) * cond1(R);
	cout << "Оценка сверху числа обусловленности (p = 1):" << markup1 << endl;*/

	/*ddouble_t markupInf = condInf(Q) * condInf(R);
	//cout << "Оценка сверху числа обусловленности (p = 1):" << markup1 << endl;*/

	////решение с возмущенной частью (deltab тут не нужна)
	//auto deltab = outrage(vec_B_outrage, 0);//внесли возмущенность 0.01(1 - "+0.01", 0 - "-0.01")
	//cout << endl << "Решение методом Гаусса с возмущенной правой частью:";
	//print_vec(gauss(vec_A_outrage, vec_B_outrage));




	//оценка снизу
	/*ddouble_t markdown1 = mark_down1(vec_A_outrage, vec_B_outrage, vec_X1);
	cout << "Оценка снизу (р = 1): " << markdown1 << endl;*/

	ddouble_t markdownInf = mark_downInf(vec_A_outrage, vec_B_outrage, vec_X1);
	cout << "Оценка снизу (р = Inf): " << markdownInf << endl;

}























