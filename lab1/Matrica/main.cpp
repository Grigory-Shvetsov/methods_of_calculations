#include<iostream>
#include<iomanip>
#include<vector>
#include<string>
#include<fstream>
#include<algorithm>

using namespace std;

// ��������� ����������!!!

//����� ������� �
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

//����� ����������� �������
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



//����� �������
template<typename T>
void print_vec(const vector<T>& vec)
{
	cout << endl;
	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << fixed << setprecision(8) << vec[i] << "\t";
	cout << endl;
}


//���������� �� �����
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


//������� �����, ���� ���� �������
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


//�������� ���
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






//������������ ����� (��������� �����)
template<typename T>
void main_str(size_t i, vector<vector<T>>& vec_A, vector<T>& vec_B)
{
	size_t k = i;
	// ���� ��������
	T max = abs(vec_A[i][i]);
	for (size_t j = i + 1; j < vec_B.size(); j++) {
		// TODO ����� ��������� �� ���� (�� ��� ������������ � �������� ���� ��� ����� � ��������� �������)
		if (abs(vec_A[j][i]) >= max) {
			max = abs(vec_A[j][i]);
			k = j;
		}
	}

	// ��� �������
	if (max <= numeric_limits<T>::epsilon() * 10)
		throw "����������� �������";

	// ������ ������
	swap(vec_A[k], vec_A[i]);
	swap(vec_B[k], vec_B[i]);
}




//����� ������
template<typename T>
vector<T> gauss(vector<vector<T>>& A, vector<T>& B) {
	// �������, ��� ������� A ���������� � ����� ����� �� ������, ��� � � B
	size_t N = A.size();
	vector<T> X;
	// ��� �������
	T c_ik;
	for (size_t k = 0; k < N - 1; ++k) {
		main_str(k, A, B);//������ ��������� � max-���������
		for (size_t i = k + 1; i < N; ++i) {
			//for (size_t i = k; i < N; ++i) {
			c_ik = A[i][k] / A[k][k];

			B[i] -= c_ik * B[k];
			for (size_t j = k + 1; j < N; ++j)
				A[i][j] -= c_ik * A[k][j];
		}
	}

	// ����� �� ��� ��� ��������� ���� (��, ��� ���� ������� A ����� ����� ����������� �������)
	if (abs(A[N - 1][N - 1]) <= numeric_limits<T>::epsilon() * 10)
		throw "����������� �������";

	reverse_steps(N, A, B, X); // �������� ���
	return X;
}



// ��������� ������� � �������
template<typename T>
vector<T> mult(const vector<vector<T>>& A1, const vector<T>& A2) {
	size_t n = A2.size();
	vector<T> res;
	res.reserve(n); // ������������ ������
	for (size_t i = 0; i < n; ++i) {
		T c = 0;
		for (size_t j = 0; j < n; ++j) {
			c += A1[i][j] * A2[j];
		}
		res.push_back(c);

	}
	return res;
}



//������ �����
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


//�������
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


// ��������� ������
template<typename T>
vector<vector<T>> mult(vector<vector<T>>& A1, vector<vector<T>>& A2) {
	size_t l = A1.size();
	size_t m = A2.size();
	size_t n = A2[0].size();

	vector<vector<T>> res;
	res.reserve(l);
	for (size_t i = 0; i < l; ++i) {
		vector<T> row;
		row.reserve(n); // ������������ ������
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


// ������� ��������� ������� �������� NxN
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



// ��������� �� ������� �������� T_ij �����
// initIndex ��������� �� ��, � ������ ������� ��������� �������� (����� ��� ������� R � QR-����������)
template<typename T>
void multT(vector<vector<T>>& A, size_t ti, size_t tj, T c, T s, size_t initIndex = 0) {
	size_t n = A.size();
	// ��������, ������ ��� � ����� ��� ������ ����� ��������
	vector<T> Ati(A[ti]);
	vector<T> Atj(A[tj]);

	for (size_t j = initIndex; j < n; ++j) {
		A[ti][j] = c * Ati[j] + s * Atj[j];
		A[tj][j] = -s * Ati[j] + c * Atj[j];
	}
}


// ���������������� ���������� ������� (�������� ������� ������� �� ����������������)
template<typename T>
void transpose(vector<vector<T>>& A) {
	size_t n = A.size();
	// ��� ���� ������
	for (size_t i = 0; i < n - 1; ++i)
		for (size_t j = i + 1; j < n; ++j)
			swap(A[i][j], A[j][i]);
}



// ������������ ������� A � QR-����������, ��� Q - ������������� �������, R - ������������������ (������������ � A)
template<typename T>
void qr(vector<vector<T>>& A, vector<vector<T>>& Q)
{
	// ��� ���� ������� ������� A ������ ���� ����������
	size_t n = A.size();
	Q = eig<T>(n);

	// ������� ������� �� ����� c � s
	for (size_t i = 0; i < n - 1; ++i)
	{

		for (size_t j = i + 1; j < n; ++j) {
			T Aii = A[i][i];
			T Aji = A[j][i];

			// ����������
			T denom = sqrt(Aii * Aii + Aji * Aji);

			if (denom <= numeric_limits<T>::epsilon() * 10)
				continue;

			T c_ij = Aii / denom;
			T s_ij = Aji / denom;

			// ��� �������
			multT(Q, i, j, c_ij, s_ij);
			multT(A, i, j, c_ij, s_ij, i);
		}
	}
	// ��� �������
	// �������� �� �������������
	for (size_t i = 0; i < n; ++i)
		if (abs(A[i][i]) <= numeric_limits<T>::epsilon() * 10)
			throw "����������� �������";
}



template<typename T>
// ��� �������
vector<T> fun_qr(vector<vector<T>>& A, vector<vector<T>>& Q, const vector<vector<T>>& vec_A_source, const vector<T>& vec_B_source)
{
	qr(A, Q);
	vector<T> vec_X2;
	reverse_steps(vec_B_source.size(), A, mult(Q, vec_B_source), vec_X2);
	/*transpose(Q);*/
	return vec_X2;
}




// ����� p=2 ���������� �������
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

// ����� p=1 ���������� �������
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

// ����� p=Inf ���������� �������
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


// ����� p=2 �������
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

// ����� p=1 �������
template <typename T>
T norm1(const vector<T>& A)
{
	size_t n = A.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		res += abs(A[i]);

	return res;
}

// ����� p=Inf �������
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


//����������
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
		throw "������������ �������� k";
	}
	return deltab;
}




// ������� �������� � ���������� �������
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
		// ������� ��������� �������
		//vector<T> e(n, 0);
		//e[i] = 1;
		// TODO ��������� �� ������� (0, 0, ..., 1, ..., 0) ����� ������� �����
		//auto B = mult(Q, e);

		// ��� �������
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



// ��������������� ������� (������ �������� �������)
template<typename T>
T cond1(vector<vector<T>>& A) {
	T norm1A = norm1Matr(A);

	auto Ainv = inverse(A);
	T norm1Ainv = norm1Matr(Ainv);

	return norm1A * norm1Ainv;
}


// ��������������� ������� (������ �������� �������)
template<typename T>
T condInf(vector<vector<T>>& A) {
	T norm1A = normInfMatr(A);

	auto Ainv = inverse(A);
	T norm1Ainv = normInfMatr(Ainv);

	return norm1A * norm1Ainv;
}


//������ ����� (� = 1)
template<typename T>
T mark_down1(const vector<vector<T>>& A, const vector<T>& B, const vector<T>& X)
{
	vector<T> deltab, deltax;
	vector<vector<T>> vec_deltab(3), vec_deltax(3);

	//��������� ������-� � ������-�
	for (size_t i = 0; i < 3; i++) {
		vector<T> B1(B);
		vec_deltab[i] = outrage(B1, i);
		//vec_deltax[i] = inverse(A) * vec_deltab[i];
		vec_deltax[i] = mult(inverse(A), vec_deltab[i]);
	}
	//��������� ��������� ������-� � ������-�
	for (size_t i = 0; i < 3; i++) {
		deltab.push_back(norm1(vec_deltab[i]) / norm1(B));
		deltax.push_back(norm1(vec_deltax[i]) / norm1(X));
	}

	vector<T> vec_deltax_deltab(3);
	//��� ������� ��������� ������-x � ������-b
	for (size_t i = 0; i < 3; i++) {
		vec_deltax_deltab[i] = deltax[i] / deltab[i];
	}
	//���������� �������� �� ���� ����������(� ������ �������� ��������� ��� � ������� ��������)
	//return max(vec_deltax_deltab);
	return *max_element(vec_deltax_deltab.begin(), vec_deltax_deltab.end());
}



//������ ����� (� = inf)
template<typename T>
T mark_downInf(const vector<vector<T>>& A, const vector<T>& B, const vector<T>& X)
{
	vector<T> deltab, deltax;
	vector<vector<T>> vec_deltab(3), vec_deltax(3);

	//��������� ������-� � ������-�
	for (size_t i = 0; i < 3; i++) {
		vector<T> B1(B);
		vec_deltab[i] = outrage(B1, i);
		//vec_deltax[i] = inverse(A) * vec_deltab[i];
		vec_deltax[i] = mult(inverse(A), vec_deltab[i]);
	}
	//��������� ��������� ������-� � ������-�
	for (size_t i = 0; i < 3; i++) {
		deltab.push_back(normInf(vec_deltab[i]) / normInf(B));
		deltax.push_back(normInf(vec_deltax[i]) / normInf(X));
	}

	vector<T> vec_deltax_deltab(3);
	//��� ������� ��������� ������-x � ������-b
	for (size_t i = 0; i < 3; i++) {
		vec_deltax_deltab[i] = deltax[i] / deltab[i];
	}
	//���������� �������� �� ���� ����������(� ������ �������� ��������� ��� � ������� ��������)
	//return max(vec_deltax_deltab);
	return *max_element(vec_deltax_deltab.begin(), vec_deltax_deltab.end());
}



// TODO �������:
// 1. ������ ����� ���������������

// �������� �������:
// 1. ������� QR-����������, ����� �� ��� ��������
// 2. ������ n ����
// 3. ��������� ��

int main()
{
	setlocale(LC_ALL, "Russian");

	typedef double ddouble_t;


	//������ �������
	vector<vector<ddouble_t>> vec_A;// ������ � �������� ������
	vector<ddouble_t> vec_B;//������ � �������� ������

	//��������� �� ����� ����������� � �������
	//size_t N = out_of_file(vec_A, vec_B, "������� 1.txt");
	size_t N = out_of_file(vec_A, vec_B, "������� 1.txt");

	// �������� ������� � ������, ������ ��� ����� �� ��������
	vector<ddouble_t> vec_B_source(vec_B);
	vector<vector<ddouble_t>> vec_A_source(vec_A);

	//����� ��� ������� fun_outrage - ������� ������� � �����������
	vector<double_t> vec_B_outrage(vec_B);
	vector<vector<ddouble_t>> vec_A_outrage(vec_A);


	cout << "����������� �������:";
	print_AB(vec_A, vec_B);


	//����� ������ (������ ������� � ������ �����)
	auto vec_X1 = gauss(vec_A, vec_B);
	cout << endl << "������� ������� ������:";
	print_vec(vec_X1); //�������


	//������� p=1
	cout << endl << "������� �=1: ";
	cout << (norm1(nevazka(vec_A_source, vec_B_source, vec_X1)));

	//������� p=2
	cout << endl << "������� �=2: ";
	cout << (norm(nevazka(vec_A_source, vec_B_source, vec_X1)));

	//������� p=inf
	cout << endl << "������� �=inf: ";
	cout << (normInf(nevazka(vec_A_source, vec_B_source, vec_X1)));



	//����� QR-���������� � ����� ������ Q � R
	vector<vector<ddouble_t>> Q;
	vector<vector<ddouble_t>> R = vec_A_source;
	auto vec_X2 = fun_qr(R, Q, vec_A_source, vec_B_source);
	cout << endl << "������� QR-�������:";
	print_vec(vec_X2); //�������

	//������� Q � R
	cout << endl << "������� Q:";
	print_A(Q);
	cout << endl << "������� R:";
	print_A(R);

	//������� p=1
	cout << endl << "������� �=1: ";
	cout << (norm1(nevazka(vec_A_source, vec_B_source, vec_X2)));

	//������� p=2
	cout << endl << "������� �=2: ";
	cout << (norm(nevazka(vec_A_source, vec_B_source, vec_X2)));

	//������� p=inf
	cout << endl << "������� �=inf: ";
	cout << (normInf(nevazka(vec_A_source, vec_B_source, vec_X2)));


	cout << endl;



	//�������� �������
	/*vector<vector<ddouble_t>> A_inv(vec_A_source);
	vector<vector<ddouble_t>> vec_A_inv = inverse(A_inv);
	cout << "�������� �������:";
	print_A(vec_A_inv);*/

	//�������� �������� �������
	/*print_A(vec_A_source);
	print_A(mult(vec_A_inv,vec_A_source));*/


	ddouble_t cond_1 = cond1(vec_A_source);
	cout << "��������������� (p = 1): " << cond_1 << endl;

	//ddouble_t cond_1 = cond1(eig<double>(4));
	//cout << "��������������� (): " << cond_1 << endl;

	/*ddouble_t cond_Inf = condInf(vec_A_source);
	cout << "��������������� (p = inf): " << cond_Inf << endl;*/



	//������ ������
	//���� ����� ���������� ������� A=A1*A2
	//������ ������ cond A <= cond A1 * cond A2
	/*ddouble_t markup1 = cond1(Q) * cond1(R);
	cout << "������ ������ ����� ��������������� (p = 1):" << markup1 << endl;*/

	/*ddouble_t markupInf = condInf(Q) * condInf(R);
	//cout << "������ ������ ����� ��������������� (p = 1):" << markup1 << endl;*/

	////������� � ����������� ������ (deltab ��� �� �����)
	//auto deltab = outrage(vec_B_outrage, 0);//������ ������������� 0.01(1 - "+0.01", 0 - "-0.01")
	//cout << endl << "������� ������� ������ � ����������� ������ ������:";
	//print_vec(gauss(vec_A_outrage, vec_B_outrage));




	//������ �����
	/*ddouble_t markdown1 = mark_down1(vec_A_outrage, vec_B_outrage, vec_X1);
	cout << "������ ����� (� = 1): " << markdown1 << endl;*/

	ddouble_t markdownInf = mark_downInf(vec_A_outrage, vec_B_outrage, vec_X1);
	cout << "������ ����� (� = Inf): " << markdownInf << endl;

}























