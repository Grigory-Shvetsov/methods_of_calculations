#include <iostream>
#include <vector>
#include <tuple>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

// TODO:
// (1). построение равномерной и чебышевской функции для заданной функции
// (2). интерполяция полиномами Лагранжа
// 3.1 метод прогонки
// 3.2 интерполяция сплайнами
// 4. построение сетки для интерполированной функции (для графиков)
// 5. сохранение этой сетки в файл
// 6. число узлов при интерполяции Лагранжа для достижения точности


template <typename T>
double mzs(T left, T right)
{
	double eps = 1e-10;
	double tao = (1 + sqrt(5)) / 2;
	double x1 = left + ((right - left) - (right - left) / tao);
	double x2 = left + (right - left) / tao;
	double foo1 = exp(x1);
	double foo2 = exp(x2);
	while (right - left > eps)
	{
		if (foo1 < foo2)
		{
			left = x1;
			x1 = x2;
			x2 = left + right - x1;
			foo1 = foo2;
			foo2 = exp(x2);
		}
		else
		{
			right = x2;
			x2 = x1;
			x1 = left + right - x2;
			foo2 = foo1;
			foo1 = exp(x1);
		}
	}
	return exp((right + left) / 2);
}

// базовый класс точки
template<typename T>
struct Point {
	T x, y;
};

template<typename T>
void print(vector<Point<T>> vec)
{
	for (const auto& p : vec) {
		cout << p.x << "\t" << p.y << endl;
	}
}


// равномерная сетка вычисляет иксы и игроки
// n - количество интервалов (будет n+1 точка)
template<typename T, typename F>
vector<Point<T>> linspace2d(F fun, T a, T b, size_t n) {
	vector<Point<T>> res;
	res.reserve(n + 1);

	T h = (b - a) / n;
	for (size_t i = 0; i <= n; ++i) {
		T x = a + h * i;
		T y = fun(x);
		res.push_back(Point<T>{ x, y });
		//res.emplace_back(x, y);
	}
	return res;
}


// равномерная сетка вычисляет только иксы
// n - количество интервалов (будет n+1 точка)
template<typename T>
vector<T> linspace1d(T a, T b, size_t n) {
	vector<T> res;
	res.reserve(n + 1);

	T h = (b - a) / n;
	for (size_t i = 0; i <= n; ++i) {
		T x = a + h * i;
		res.push_back(x);
	}
	return res;
}


// чебышевская сетка вычисляет иксы и игроки
template<typename T, typename F>
vector<Point<T>> chebyshevGrid2d(F fun, T a, T b, size_t n) {
	vector<Point<T>> res;
	res.reserve(n + 1);

	// проходим цикл в обратном порядке, потому что формула сетки их так выдает
	// защита от переполнения
	for (size_t i = n; i >= 0 && i <= n; --i) {
		T x = (a + b) / 2.0 + (b - a) / 2.0 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
		T y = fun(x);
		res.push_back(Point<T>{ x, y });
	}
	return res;
}


// чебышевская сетка вычисляет только иксы
template<typename T>
vector<T> chebyshevGrid1d(T a, T b, size_t n) {
	vector<T> res;
	res.reserve(n + 1);

	// проходим цикл в обратном порядке, потому что формула сетки их так выдает
	// защита от переполнения
	for (size_t i = n; i >= 0 && i <= n; --i) {
		T x = (a + b) / 2.0 + (b - a) / 2.0 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
		res.push_back(x);
	}
	return res;
}


// интерполяция полиномами Лагранжа
template<typename T>
class LagrangePolynomial {
	vector<Point<T>> mPoints;
	size_t n;
public:
	LagrangePolynomial(const vector<Point<T>>& points) : mPoints(points) {
		n = points.size() - 1;
	}
	T operator()(T x) {
		T sum = 0;
		for (size_t k = 0; k <= n; ++k) {
			T ck = 1;
			T xk = mPoints[k].x;
			//обходим нулевой знаменатель
			for (size_t j = 0; j < k; ++j) {
				T xj = mPoints[j].x;
				ck *= (x - xj) / (xk - xj);
			}
			for (size_t j = k + 1; j <= n; ++j) {
				T xj = mPoints[j].x;
				ck *= (x - xj) / (xk - xj);
			}
			sum += mPoints[k].y * ck;
		}
		return sum;
	}


	// TODO можно сделать доп. функции вроде оценки погрешности, максимальный шаг сетки и т.д.
};

// TODO метод прогонки
template<typename T>
vector<T> ThomasAlgorithm(
	const vector<T>& a, const vector<T>& b, const vector<T>& c,
	const vector<T>& d
) {
	size_t n = b.size();
	vector<T> alpha, beta;
	alpha.reserve(n - 1);
	beta.reserve(n - 1);

	// i = 1
	alpha.push_back(-c[0] / b[0]);
	beta.push_back(d[0] / b[0]);
	// i = 2..n-1
	for (size_t i = 1; i < n - 1; ++i) {
		T denom = -b[i] - a[i - 1] * alpha.back();
		alpha.push_back(c[i] / denom);
		beta.push_back((-d[i] + a[i - 1] * beta.back()) / denom);
	}
	vector<T> X;
	X.reserve(n);
	X.push_back((-d[n - 1] + a[n - 2] * beta[n - 2]) / (-b[n - 1] - a[n - 2] * alpha[n - 2]));
	for (size_t i = n - 2; i >= 0 && i < n; --i)
		X.push_back(alpha[i] * X.back() + beta[i]);

	reverse(X.begin(), X.end());
	return X;
}


// интерполяция кубическими сплайнами
template<typename T>
class SplineInterpolation {
	// сохраняем иксы для определения нужно интервала
	vector<T> pX;
	vector<T> a, b, c, d;
	size_t n;

	// ищет номер необходимого интервала с помощью бинарного поиска
	// Если что-то пошло не так (что?), возвращает SIZE_MAX
	size_t binarySearch(T x) {
		// Если лежим вне отрезка, то возвращаем ближайший (первый или последний)
		if (x < pX[0])
			return 0;
		else if (pX[n] < x)
			return n - 1;

		size_t low = 0;
		//size_t high = n; // TODO не должно быть n-1?
		size_t high = n - 1;
		while (low <= high) {
			size_t mid = (low + high) / 2;
			T xm = pX[mid];
			T xm1 = pX[mid + 1];
			if (xm1 < x)
				low = mid + 1;
			else if (x < xm)
				high = mid - 1;
			else
				return mid;
		}
		return SIZE_MAX;
	}
public:
	SplineInterpolation(const vector<Point<T>>& points) : n(points.size() - 1) {
		pX.reserve(points.size());
		a.reserve(n);
		for (size_t i = 0; i < n; ++i) {
			Point<T> p = points[i];
			pX.push_back(p.x);
			a.push_back(p.y);
		}
		pX.push_back(points.back().x);

		vector<T> h, g;
		h.reserve(n);
		g.reserve(n);
		for (size_t i = 1; i <= n; ++i) {
			Point<T> pi1 = points[i - 1];
			Point<T> pi = points[i];
			T hi = pi.x - pi1.x;
			T gi = (pi.y - pi1.y) / hi;
			h.push_back(hi);
			g.push_back(gi);
		}
		// вычисляем коэффиценты
		// составляем СЛАУ
		vector<T> aS(h.begin() + 1, h.end() - 1);
		vector<T> bS;
		bS.reserve(n - 1);
		vector<T> cS(h.begin() + 1, h.end() - 1);
		vector<T> dS;
		dS.reserve(n - 1);
		for (size_t i = 1; i < n; ++i) {
			bS.push_back(2 * (h[i - 1] + h[i]));
			dS.push_back(3 * (g[i] - g[i - 1]));
		}
		// решаем
		c = ThomasAlgorithm(aS, bS, cS, dS);
		c.insert(c.begin(), 0);
		//c.push_back(0);

		// заполняем остальные коэффициенты
		b.reserve(n);
		d.reserve(n);
		for (size_t i = 0; i < n - 1; ++i) {
			b.push_back(g[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3);
			d.push_back((c[i + 1] - c[i]) / (3 * h[i]));
		}

		b.push_back(g[n - 1] - 2 * c[n - 1] * h[n - 1] / 3);
		d.push_back(-c[n - 1] / (3 * h[n - 1]));
	}

	// TODO можно сделать бинарный поиск нужно интервала поиск
	T operator()(T x) {
		/*for (size_t i = 0; i < n; ++i)
			if (pX[i] <= x && x <= pX[i + 1]) {
				T xi = pX[i];
				T xxi = x - xi;
				// (x-xi)^2
				T xxi2 = xxi * xxi;
				// (x-xi)^3
				T xxi3 = xxi2 * xxi;
				return a[i] + b[i] * xxi + c[i] * xxi2 + d[i] * xxi3;
			}
		throw "x не принадлежит диапозону интерполяции";*/

		size_t i = binarySearch(x);
		//cout << "i == " << i << endl;
		if (i == SIZE_MAX)
			throw "x не принадлежит диапозону интерполяции";

		T xi = pX[i];
		T xxi = x - xi;
		// (x-xi)^2
		T xxi2 = xxi * xxi;
		// (x-xi)^3
		T xxi3 = xxi2 * xxi;
		return a[i] + b[i] * xxi + c[i] * xxi2 + d[i] * xxi3;
	}
};


template<typename T>
ostream& operator<<(ostream& out, const vector<Point<T>>& points) {
	for (auto& point : points)
		out << point.x << "\t" << point.y << endl;
	return out;
}

template<typename T>
T x2(T x) {
	return x * x;
}

template<typename T>
T x4(T x) {
	return x * x * x * x;
}

template<typename T>
T xx(T x) {
	return x;
}

template<typename T>
T msin(T x) {
	return sin(x);
}

/////// Примеры из методички

// -1 <= x <= 1
template<typename T>
T example1(T x) {
	return x * x;
}

// пример Рунге, -1 <= x <= 1
template<typename T>
T example2(T x) {
	return 1.0 / (1 + x * x);
}

// -3 <= x <= 3
template<typename T>
T example3(T x) {
	return 1.0 / atan(1.0 + 10.0 * x * x);
}

// -1 <= x <= 1
template<typename T>
T example4(T x) {
	return pow(4.0 * x * x * x + 2.0 * x * x - 4.0 * x + 2.0, sqrt(2.0)) + asin(1.0 / (5.0 + x - x * x)) - 5.0;
}


// -1 <= x <= 1
template<typename T>
T example5(T x) {
	return 1.0 / (1 + 25.0 * x * x);
}

template<typename T>
T myexp(T x) {
	return exp(x);
}


template<typename T>
T example6(T x) {
	return 1;
}



template<typename T>
T example7(T x) {
	return abs(x);
}

template <typename T>
pair<T, T> task3(T x, T a, T b, size_t n)
{
	vector<Point<T>> ps;
	ps = linspace2d(myexp<T>, a, b, n);
	LagrangePolynomial<T> task(ps);
	T mistake = abs(task(x) - exp(x));
	T grade = exp((n + 1) * log((b - a) / n)) * mzs(a, x);
	return make_pair(mistake, grade);
}


template <typename T, typename F>
tuple <vector<Point<T>>, vector<Point<T>>, vector<Point<T>>, vector<Point<T>>> calculate(
	F fun,
	const vector<Point<T>>& ps1,
	const vector<T>& ps1Plot,
	const vector<Point<T>>& ps2,
	const vector<T>& ps2Plot)
{
	//индекс 1 - равномерная сетка, 2 - чебышевская
	vector<Point<T>> resL1, resL2, resCh1, resCh2;
	LagrangePolynomial<T> testLagrange1(ps1);
	LagrangePolynomial<T> testLagrange2(ps2);
	SplineInterpolation<T> testSpline1(ps1);
	SplineInterpolation<T> testSpline2(ps2);
	Point<T> point;
	T sum1, sum2, sum3, sum4;
	sum1 = sum2 = sum3 = sum4 = 0;
	for (size_t i = 0; i < ps1Plot.size(); i++) {
		point.x = ps1Plot[i];
		point.y = testLagrange1(point.x);
		resL1.push_back(point);
		if (sum1 < abs(point.y - fun(point.x)))
			sum1 = abs(point.y - fun(point.x));

		point.x = ps2Plot[i];
		point.y = testLagrange2(point.x);
		resL2.push_back(point);
		if (sum2 < abs(point.y - fun(point.x)))
			sum2 = abs(point.y - fun(point.x));

		point.x = ps1Plot[i];
		point.y = testSpline1(point.x);
		resCh1.push_back(point);
		if (sum3 < abs(point.y - fun(point.x)))
			sum3 = abs(point.y - fun(point.x));

		point.x = ps2Plot[i];
		point.y = testSpline2(point.x);
		resCh2.push_back(point);
		if (sum4 < abs(point.y - fun(point.x)))
			sum4 = abs(point.y - fun(point.x));

	}
	cout << fixed << setprecision(15) << "Норма ошибки Лагранжа для равномерной: " << sum1 << endl;
	cout << fixed << setprecision(15) << "Норма ошибки Лагранжа для чебышевской: " << sum2 << endl;
	cout << fixed << setprecision(15) << "Норма ошибки Сплайн для равномерной: " << sum3 << endl;
	cout << fixed << setprecision(15) << "Норма ошибки Сплайн для чебышевской: " << sum4 << endl;


	return make_tuple(resL1, resL2, resCh1, resCh2);
}

template <typename T, typename F>
tuple<vector<Point<T>>, vector<T>, vector<Point<T>>, vector<T>> grid(
	F fun,
	T a,
	T b,
	size_t n,
	size_t nPlot)
{
	vector<Point<T>> ps1, ps2;
	vector<T> ps1Plot, ps2Plot;
	//создаем две сетки и два вектора иксов
	ps1 = linspace2d(fun, a, b, n);
	ps1Plot = linspace1d(a, b, nPlot);
	ps2 = chebyshevGrid2d(fun, a, b, n);
	ps2Plot = chebyshevGrid1d(a, b, nPlot);
	return make_tuple(ps1, ps1Plot, ps2, ps2Plot);
}

template <typename T>
void inFile(
	const vector<Point<T>>& resL1,
	const vector<Point<T>>& resL2,
	const vector<Point<T>>& resCh1,
	const vector<Point<T>>& resCh2)
{
	ofstream outL1("Л_Рсетка.txt");
	ofstream outL2("Л_Чсетка.txt");
	ofstream outCh1("С_Рсетка.txt");
	ofstream outCh2("С_Чсетка.txt");
	for (int i = 0; i < resL1.size(); i++) {
		outL1 << resL1[i].x << "\t" << resL1[i].y << endl;
	}
	for (int i = 0; i < resL2.size(); i++) {
		outL2 << resL2[i].x << "\t" << resL2[i].y << endl;
	}
	for (int i = 0; i < resCh1.size(); i++) {
		outCh1 << resCh1[i].x << "\t" << resCh1[i].y << endl;
	}
	for (int i = 0; i < resCh2.size(); i++) {
		outCh2 << resCh2[i].x << "\t" << resCh2[i].y << endl;
	}
	outL1.close();
	outL2.close();
	outCh1.close();
	outCh2.close();

}

int main()
{
	typedef double MyType;
	setlocale(LC_ALL, "Russian");
	size_t n = 3;
	size_t nPlot = 10000;
	vector<Point<MyType>> ps1, ps2;
	vector<MyType> ps1Plot, ps2Plot;
	//меняем функцию просто в grid
	tie(ps1, ps1Plot, ps2, ps2Plot) = grid(x2<MyType>, -1.0, 1.0, n, nPlot);


	vector<Point<MyType>> resL1, resL2, resCh1, resCh2;
	tie(resL1, resL2, resCh1, resCh2) = calculate(x2<MyType>, ps1, ps1Plot, ps2, ps2Plot);
	//print(res1);
	cout << endl;
	//print(res2);
	inFile(resL1, resL2, resCh1, resCh2);

	//Для 3 контрольного вопроса
	/*MyType mistake, grade;
	tie(mistake, grade) = task3(2.2, 0.0, 2.0, 10);
	cout << "Ошибка |exp - L_n| = " << mistake << " <= " << grade << endl;;*/


}