#include "function.cpp"
#include "helpfun.cpp"


//double t7exakt = 0.516127723925248;
double t7exakt = 1.0;

// TODO:
// (1). Локализация корня (нужно ли это делать для системы? и как это тогда делать?)
// (2). Метод бисекции
// (3). Метод Ньютона
// 4. Модификация метода Ньютона (с методом хорд)
// 5. Модификация метода Ньютона для системы уравнений (передавать матрицу Якоби или вычислять её численно)
// 5*. Просто метод Ньютона для системы (для диаграммы сходимости)
// (5**). Обращение матрицы
// (6). Численное вычисление матрицы Якоби

// Локализация корня
// Принимает:
// fun - функция, которой происходит локализация корня
// a, b - отрезок локализации
// n - кол-во подотрезков (по умолчанию 100)
// Возвращает: вектор отрезков, которые содержат корни
template<typename T, typename F>
vector<pair<T, T>> localize(F fun, T a, T b, size_t n = 100) {
	vector<pair<T, T>> res;
	T h = (b - a) / n;
	T f0 = fun(a), f1;
	T x0 = a, x1;
	for (size_t i = 1; i <= n; ++i) {
		x1 = a + h * i;
		f1 = fun(x1);

		//cout << x0 << " " << x1 << endl;

		if (f0 * f1 < 0)
			res.push_back(make_pair(x0, x1));

		//cout << x0 << "->" << f0 << "; " << x1 << "->" << f1 << endl;

		x0 = x1;
		f0 = f1;
	}
	return res;
}

// локализация корней
template <typename T, typename F>
vector<pair<T, T>> localization(F fun, T a, T b, size_t n)
{
	T h = (b - a) / n;
	vector<T> setka;
	for (size_t i = 0; i < n; i++) {
		setka.push_back(a + h * i);
	}
	setka.push_back(b);
	//print_vec(setka);
	pair<T, T> local;
	vector<pair<T, T>> locals;
	vector<T> roots;
	size_t i = 0;
	while (i < n)
	{
		if (fun(setka[i]) * fun(setka[i + 1]) < 0.0)
		{
			local.first = setka[i];
			local.second = setka[i + 1];
			locals.push_back(local);
		}
		else if (fun(setka[i]) * fun(setka[i + 1]) == 0.0) // если значений из сетки и есть корень
		{
			if (fun(setka[i]) == 0.0)
				roots.push_back(setka[i]);
			else
				roots.push_back(setka[i + 1]);
		}
		++i;
	}

	return locals;
}

// Метод бисекции
template<typename T, typename F>
T bisection(F fun, T a, T b, T eps, size_t& iterCount) {
	T ak = a;
	T bk = b;

	T fak = fun(ak);
	T fbk = fun(bk);
	T xk = (a + b) / 2;
	vector<T> X;
	X.push_back(xk);
	while (abs(bk - ak) > 2 * eps) {
		T fxk = fun(xk);
		if (fak * fxk <= 0) {
			bk = xk;
			fbk = fxk;
		}
		else {
			ak = xk;
			fak = fxk;
		}

		xk = (ak + bk) / 2;
		X.push_back(xk);
		++iterCount;
	}
	degreeX(X, t7exakt);
	return xk;
}

template<typename T, typename F>
vector<T> bisection(F fun, vector<pair<T, T>> intervals, T eps) {
	vector<T> res;
	res.reserve(intervals.size());
	size_t iterCount = 0;
	for (pair<T, T>& interval : intervals) {
		T a = interval.first;
		T b = interval.second;
		T x = bisection(fun, a, b, eps, iterCount);
		res.push_back(x);
	}
	cout << "Количество итераций для нахождения всех корней:" << iterCount << endl;
	return res;
}

// Метод Ньютона (без модификации)
// diff - производная функции fun
template<typename T, typename F, typename Diff>
T Newton(F fun, Diff diff, T a, T b, T x0, T eps, size_t& iterCount, size_t p = 1) {
	T xk = x0;
	T coef;
	vector<T> X;
	X.push_back(xk);
	do {
		T newXk = xk - p * fun(xk) / diff(xk);
		coef = abs(newXk - xk);
		xk = newXk;
		X.push_back(xk);
		++iterCount;
		//cout << newXk <<"\t" << abs(newXk - t7exakt) << endl;
	} while (coef > eps);
	degreeX(X, t7exakt);
	return xk;
}


// Упрощенный метод Ньютона
template<typename T, typename F, typename Diff>
T NewtonX0(F fun, Diff diff, T a, T b, T x0, T eps, size_t& iterCount) {
	T xk = x0;
	T coef;
	T diffx0 = diff(x0);
	do {
		T newXk = xk - fun(xk) / diffx0;
		coef = abs(newXk - xk);
		xk = newXk;
		++iterCount;
	} while (coef > eps);
	return xk;
}



// Модификация метода Ньютона с защитой от выхода за границу
template<typename T, typename F, typename Diff>
T NewtonMod(F fun, Diff diff, T a, T b, T x0, T eps) {
	T fa = fun(a);
	T fb = fun(b);

	T xk = x0;
	T coef;

	size_t iterCount = 0;

	do {
		T newXk;
		if ((iterCount + 1) % 9 != 0)
			newXk = xk - fun(xk) / diff(xk);
		else {
			// Защита от зацикливания (просто раз в 9 итераций применяем
			// другой метод, в данном случае метод бисекции)
			T fxk = fun(xk);


			// метод бисекции
			if (fa * fxk <= 0)
				newXk = (a + xk) / 2;
			else
				newXk = (xk + b) / 2;

			// модифицированный метод хорд
			/*if (fa * fxk <= 0)
				newXk = xk - fxk * (a - xk) / (fa - fxk);
			else
				newXk = xk - fxk * (b - xk) / (fb - fxk);*/
		}

		// Защита модифицированным методом бисекции
		if (newXk < a || newXk > b) {
			T fxk = fun(xk);

			// метод бисекции
			/*if (fa * fxk <= 0)
				newXk = (a + xk) / 2;
			else
				newXk = (xk + b) / 2;
				*/

				// модифицированный метод хорд
			if (fa * fxk <= 0)
				newXk = xk - fxk * (a - xk) / (fa - fxk);
			else
				newXk = xk - fxk * (b - xk) / (fb - fxk);

		}

		coef = abs(newXk - xk);
		xk = newXk;
		iterCount++;
	} while (coef > eps);

	cout << "Количество итераций: " << iterCount << endl;

	return xk;
}




// Метод Ньютона (без модификации)
// diff - производная функции fun
template<typename T, typename F, typename Diff>
T Newton(F fun, Diff diff, T a, T b, T eps, size_t& iterCount, int method) {
	// Начальная итерация методом хорд
	T fa = fun(a);
	T fb = fun(b);
	T x0 = (fa * b - fb * a) / (fa - fb);
	++iterCount;
	switch (method)
	{
	case 0:		return Newton(fun, diff, a, b, x0, eps, iterCount);//классика
	case 1: 	return NewtonX0(fun, diff, a, b, x0, eps, iterCount);//f'(x0)
	case 2: 	return NewtonMod(fun, diff, a, b, x0, eps);// модификация
	default:
		break;
	}
}


template<typename T, typename F, typename Diff>
vector<T> Newton(F fun, Diff diff, vector<pair<T, T>> intervals, T eps, int method) {
	vector<T> res;
	res.reserve(intervals.size());
	size_t iterCount = 0;
	for (pair<T, T>& interval : intervals) {
		T a = interval.first;
		T b = interval.second;
		T x = Newton(fun, diff, a, b, eps, iterCount, method);
		res.push_back(x);
	}
	cout << "Количество итераций для нахождения всех корней: " << iterCount << endl;
	return res;
}





// Классический метод Ньютона для системы уравенений
// На каждой итерации считается матрица Якоби
template<typename T, typename Fun, typename J>
//vector<T> NewtonSystem(Fun f1, Fun f2, J jacobi, T l1, T l2, T eps, const vector<T>& x0)
vector<T> NewtonSystem(Fun f1, Fun f2, J jacobi, T eps, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

	size_t iterCount = 0;

	do {
		// F'^-1
		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));
		vector<T> F{ f1(xk[0], xk[1]), f2(xk[0], xk[1]) };

		// F'^-1 * F
		vector<T> jacobiInvMultF = mult(jacobiInv, F);

		// x^k+1 = x^k - F'^-1 * F
		// x^k+1 - x^k = - F'^-1 * F
		// |x^k+1 - x^k| = |- F'^-1 * F| = |F'^-1 * F|
		// => coef = norm(jacobiInvMultF)
		vector<T> newXk;
		// Типо задел под n-мерный случай
		newXk.reserve(x0.size());
		for (size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		iterCount++;

		/*if (abs(xk[0]) > l1 || abs(xk[1]) > l2) {
			cout << "Вышли за границу ";
			print_vec(xk);
		}*/

	} while (coef > eps);

	cout << "Количество итераций: " << iterCount << endl;

	return xk;
}





// Классический метод Ньютона для системы уравенений
// Матрица Якоби считается только для первого приближения
template<typename T, typename Fun, typename J>
//vector<T> NewtonSystem(Fun f1, Fun f2, J jacobi, T l1, T l2, T eps, const vector<T>& x0)
vector<T> NewtonSystemX0(Fun f1, Fun f2, J jacobi, T eps, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;
	matrix<T> jacobiInv = inv(jacobi(x0[0], x0[1]));
	size_t iterCount = 0;

	do {
		// F'^-1
		vector<T> F{ f1(xk[0], xk[1]), f2(xk[0], xk[1]) };

		// F'^-1 * F
		vector<T> jacobiInvMultF = mult(jacobiInv, F);

		// x^k+1 = x^k - F'^-1 * F
		// x^k+1 - x^k = - F'^-1 * F
		// |x^k+1 - x^k| = |- F'^-1 * F| = |F'^-1 * F|
		// => coef = norm(jacobiInvMultF)
		vector<T> newXk;
		// Типо задел под n-мерный случай
		newXk.reserve(x0.size());
		for (size_t i = 0; i < x0.size(); ++i) {
			newXk.push_back(xk[i] - jacobiInvMultF[i]);
		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		iterCount++;

		/*if (abs(xk[0]) > l1 || abs(xk[1]) > l2) {
			cout << "Вышли за границу ";
			print_vec(xk);
		}*/

	} while (coef > eps);

	cout << "Количество итераций: " << iterCount << endl;

	return xk;
}




//  метод хорд
template<typename T, typename F>
T chord(F f, T a, T b, T eps, size_t& iterCount)
{
	T xk = a;
	T xk1 = b;
	T xk2, coef;
	vector<T> X;
	X.push_back(xk);
	X.push_back(xk1);
	do
	{
		++iterCount;
		xk2 = xk - f(xk) * (xk1 - xk) / (f(xk1) - f(xk));
		X.push_back(xk2);
		coef = abs(xk2 - xk1);
		xk = xk1;
		xk1 = xk2;
	} while (coef > eps);
	degreeX(X, t7exakt);
	return xk2;
}

//////////////////////



// Численная матрица Якоби
template<typename T, typename F>
auto numJacobi(F f1, F f2, T eps) {
	return [&f1, &f2, eps](T x1, T x2) {
		T f1Val = f1(x1, x2);
		T f2Val = f2(x1, x2);
		return matrix<T>{
			{ (f1(x1 + eps, x2) - f1Val) / eps, (f1(x1, x2 + eps) - f1Val) / eps },
			{ (f2(x1 + eps, x2) - f2Val) / eps, (f2(x1, x2 + eps) - f2Val) / eps }
		};
	};
}


// Модификация метода Ньютона для системы уравенений с проверкой выхода за границу
// Границы задаются следующим образом: x1 <= x <= x2, y1 <= y <= y2
template<typename T, typename Fun, typename J>
vector<T> NewtonSystemMod(Fun f1, Fun f2, J jacobi, T x1, T x2, T y1, T y2, T eps, const vector<T>& x0)
{
	vector<T> xk = x0;
	T coef;

	size_t iterCount = 0;

	do {
		// F'^-1
		matrix<T> jacobiInv = inv(jacobi(xk[0], xk[1]));
		vector<T> F{ f1(xk[0], xk[1]), f2(xk[0], xk[1]) };

		// F'^-1 * F
		vector<T> jacobiInvMultF = mult(jacobiInv, F);

		vector<T> newXk = xk;
		newXk[0] -= jacobiInvMultF[0];
		newXk[1] -= jacobiInvMultF[1];

		if (newXk[0] < x1 || newXk[0] > x2 || newXk[1] < y1 || newXk[1] > y2) {
			// Для защиты от выхода за область применяем гибридный метод:
			// Внешние итерации - одна итерация нелинейного метода Зейделя
			// Внутренние - одна итерация метода бисекции

			//cout << "Вышли за границу" << endl;

			newXk = xk;

			T fax = f1(x1, newXk[1]);
			T fbx = f1(x2, newXk[1]);

			T fxk = f1(newXk[0], newXk[1]);

			if (fax * fxk <= 0)
				newXk[0] = (x1 + newXk[0]) / 2;
			else
				newXk[0] = (x2 + newXk[0]) / 2;

			T fay = f2(newXk[0], y1);
			T fby = f2(newXk[0], y2);

			T fyk = f2(newXk[0], newXk[1]);

			if (fay * fby <= 0)
				newXk[1] = (y1 + newXk[1]) / 2;
			else
				newXk[1] = (y2 + newXk[1]) / 2;

		}

		coef = norm(jacobiInvMultF);

		xk = newXk;

		iterCount++;

	} while (coef > eps);

	cout << "Количество итераций: " << iterCount << endl;

	return xk;
}


int main()
{
	typedef double Type;
	setlocale(LC_ALL, "Russian");

	vector<Type> x0test4{ -2.433, -9.654 };
	vector<Type> x0test5{ -1., 4. };
	vector<Type> x0test6{ 0.9, 0.5 };





	//локализация
	auto intervals = localization(test2<Type>, -1.0, 10.0, 15);
	cout << "Локализация корней:" << endl;
	cout << intervals << endl;
	cout << endl;


	cout << "Метод бисекции:" << endl;
	vector<Type> rootsBis = bisection(test2<Type>, intervals, 1e-3);
	//cout << "Корни:" << endl;
	for (Type& root : rootsBis)
		cout << "x = " << root << "\t";

	cout << endl;
	cout << endl;

	cout << "Классический метод Ньютона:" << endl;
	vector<Type> rootsNew = Newton(test2<Type>, numDeriv(test2<Type>, 1e-3), intervals, 1e-3, 0);
	//cout << "Корни:" << endl;
	for (Type& root : rootsNew)
		cout << "x = " << root << "\t";


	cout << endl;
	cout << endl;


	cout << "Упрощенный метод Ньютона:" << endl;
	vector<Type> rootsNewXO = Newton(test2<Type>, numDeriv(test2<Type>, 1e-3), intervals, 1e-3, 1);
	//cout << "Корни:" << endl;
	for (Type& root : rootsNewXO)
		cout << "x = " << root << "\t";

	cout << endl;
	cout << endl;


	size_t iterCount = 0;
	cout << "Тест 2 с выходом за границы:" << endl;
	// тут зацикливание
	Type rootNewTest2 = Newton(test2<Type>, test2Df<Type>, -1., 10., 8., 1e-3, iterCount);
	// Выводит nan
	cout << "x = " << rootNewTest2 << endl;
	Type rootNewTest2Mod = NewtonMod(test2<Type>, test2Df<Type>, -1., 10., 8., 1e-3);
	// Выводит правильный ответ
	cout << "x = " << rootNewTest2Mod << endl;

	cout << endl;

	iterCount = 0;
	cout << "Тест 3 с зацикливанием:" << endl;
	// тут зацикливание
	/*Type rootNewTest3 = Newton(test3<Type>, test3Df<Type>, 0., 1., 0., 1e-3, iterCount);
	cout << "x = " << rootNewTest3 << endl;	*/
	// а тут всё ок
	Type rootNewTest3Mod = NewtonMod(test3<Type>, test3Df<Type>, 0., 1., 0., 1e-3);
	cout << "x = " << rootNewTest3Mod << endl;
	cout << endl;
	cout << endl;
	cout << "---------------------------------СИСТЕМЫ---------------------------------";
	cout << endl;
	cout << endl;


	//ТЕСТ 4

	cout << "Классический метод Ньютона ТЕСТ 4" << endl;
	vector<Type> rootNewSys = NewtonSystem(test4f1<Type>, test4f2<Type>,
		//numJacobi<Type>(test4f1<Type>, test4f2<Type>, 1e-3),
		test4Jacobi<Type>,
		// Если выбрать x0 = {0, 0}, то метод разойдется
		//10., 10., 1e-3, vector<Type>{-5.0, 2.0}
		//10., 10., 1e-3, vector<Type>{1.5, 4.3}
		1e-3,
		x0test4);
	print_vec(rootNewSys);

	cout << "Упрощенный метод Ньютона ТЕСТ 4" << endl;
	vector<Type> rootNewSysX04 = NewtonSystemX0(
		test4f1<Type>,
		test4f2<Type>,
		test4Jacobi<Type>,
		1e-3,
		x0test4
	);
	print_vec(rootNewSysX04);


	cout << "Модифицированный метод Ньютона ТЕСТ 4" << endl;
	vector<Type> rootNewSysMod4 = NewtonSystemMod(test4f1<Type>, test4f2<Type>,
		//numJacobi<Type>(test4f1<Type>, test4f2<Type>, 1e-3),
		test4Jacobi<Type>,
		-10., 10., -10., 10.,
		1e-3,
		// Если выбрать x0 = {0, 0}, то метод разойдется
		x0test4
	);
	print_vec(rootNewSysMod4);



	//ТЕСТ 5

	cout << "Классический метод Ньютона ТЕСТ 5" << endl;
	vector<Type> rootNewSys5 = NewtonSystem(
		test5f1<Type>,
		test5f2<Type>,
		test5Jacobi<Type>,
		1e-3,
		x0test5
	);
	print_vec(rootNewSys5);

	cout << "Упрощенный метод Ньютона ТЕСТ 5" << endl;
	vector<Type> rootNewSysX05 = NewtonSystemX0(
		test5f1<Type>,
		test5f2<Type>,
		test5Jacobi<Type>,
		1e-3,
		x0test5
	);
	print_vec(rootNewSysX05);

	//cout << "Модифицированный метод Ньютона ТЕСТ 5" << endl;
	//vector<Type> rootNewSysMod5 = NewtonSystemMod(test5f1<Type>, test5f2<Type>,
	//	//numJacobi<Type>(test4f1<Type>, test4f2<Type>, 1e-3),
	//	test4Jacobi<Type>,
	//	-10., 10., -10., 10.,
	//	1e-3,
	//	// Если выбрать x0 = {0, 0}, то метод разойдется
	//	x0test5
	//);
	//print_vec(rootNewSysMod5);

	//ТЕСТ 6

	cout << "Классический метод Ньютона ТЕСТ 6" << endl;
	vector<Type> rootNewSys6 = NewtonSystem(
		test6f1<Type>,
		test6f2<Type>,
		test6Jacobi<Type>,
		1e-4,
		x0test6);
	print_vec(rootNewSys6);


	cout << "Упрощенный метод Ньютона ТЕСТ 6" << endl;
	vector<Type> rootNewSysX06 = NewtonSystemX0(
		test6f1<Type>,
		test6f2<Type>,
		test6Jacobi<Type>,
		1e-4,
		x0test6);
	print_vec(rootNewSysX06);





	//cout << "Модифицированный метод Ньютона ТЕСТ 6" << endl;
	//vector<Type> rootNewSysMod6 = NewtonSystemMod(test6f1<Type>, test6f2<Type>,
	//	//numJacobi<Type>(test4f1<Type>, test4f2<Type>, 1e-3),
	//	test4Jacobi<Type>,
	//	-10., 10., -10., 10., 1e-4,
	//	// Если выбрать x0 = {0, 0}, то метод разойдется
	//	//vector<Type>{-5.0, 2.0}
	//	//vector<Type>{1.5, 4.3}
	//	//vector<Type>{-2.433, -9.654}
	//	x0test6
	//);
	//print_vec(rootNewSysMod6);



	cout << endl;
	cout << endl;
	cout << endl;

	cout << "Тест для отчета:" << endl;
	cout << fixed;
	cout.precision(15);
	auto f7 = test7<Type>;
	auto f7Df = test7Df<Type>;
	Type eps = 1e-3;
	auto f7DfNum = numDeriv(f7, eps);



	cout << "Метод хорд:" << endl;
	//size_t iterCount = 0;
	iterCount = 0;
	Type t7Chord = chord(f7, -1., 1., eps, iterCount);
	cout << "Количество итераций: " << iterCount << endl;
	cout << "Результат: " << t7Chord << endl;
	cout << "Достигнутая точность: " << abs(t7exakt - t7Chord) << endl;
	cout << "Невязка: " << abs(f7(t7Chord)) << endl << endl;


	cout << "Метод бисекции:" << endl;
	//size_t iterCount = 0;
	iterCount = 0;
	Type t7Bisection = bisection(f7, -1., 1., eps, iterCount);
	cout << "Количество итераций: " << iterCount << endl;
	cout << "Результат: " << t7Bisection << endl;
	cout << "Достигнутая точность: " << abs(t7exakt - t7Bisection) << endl;
	cout << "Невязка: " << abs(f7(t7Bisection)) << endl;
	// TODO оценка трудоемкости

	cout << endl;

	cout << "Метод Ньютона, аналитическая производная:" << endl;
	iterCount = 0;
	Type t7NewAn = Newton(f7, f7Df, -1., 1., -0.7, eps, iterCount);
	cout << "Количество итераций: " << iterCount << endl;
	cout << "Результат: " << t7NewAn << endl;
	cout << "Достигнутая точность: " << abs(t7exakt - t7NewAn) << endl;
	cout << "Невязка: " << abs(f7(t7NewAn)) << endl;
	// TODO оценка трудоемкости

	cout << endl;

	cout << "Метод Ньютона, численная производная:" << endl;
	iterCount = 0;
	Type t7NewNum = Newton(f7, numDeriv(f7, eps), -1., 1., -0.7, eps, iterCount);
	cout << "Количество итераций: " << iterCount << endl;
	cout << "Результат: " << t7NewNum << endl;
	cout << "Достигнутая точность: " << abs(t7exakt - t7NewNum) << endl;
	cout << "Невязка: " << abs(f7(t7NewNum)) << endl;
	// TODO оценка трудоемкости

	cout << "------------------------------------" << endl;

	cout << "Метод Ньютона с корнем кратности 2:" << endl;
	iterCount = 0;
	Type t0NewAn = Newton(test0<Type>, test0Df<Type>, 0., 2., 0.2, eps, iterCount, 2);
	cout << "Количество итераций: " << iterCount << endl;
	cout << "Результат: " << t0NewAn << endl;
	cout << "Достигнутая точность: " << abs(t7exakt - t0NewAn) << endl;
	cout << "Невязка: " << abs(f7(t0NewAn)) << endl;


}


