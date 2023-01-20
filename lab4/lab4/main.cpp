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

//������� �����, ���� ���� �������
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

//�������� ���
// ��������� O(n^2 / 2)
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

	// ��� ����������� ����������� � ������������ ����� ������� ����������
	// ������� �� �����������. ������� �������� ����� ����������������
	// 
	// ����� �� ��� ��� ��������� ���� (��, ��� ���� ������� A ����� ����� ����������� �������)
	//if (abs(A[N - 1][N - 1]) <= numeric_limits<T>::epsilon() * 10)
	//	throw "����������� �������";

	reverse_steps(N, A, B, X); // �������� ���
	return X;
}


// ��������� ������� � �������
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

// ��������� ������������ ���� ��������
template<typename T>
T mult(const vector<T>& A1, const vector<T>& A2) {
	size_t n = A1.size();
	T res = 0;
	for (size_t i = 0; i < n; ++i)
		res += A1[i] * A2[i];
	return res;
}


// ��������� ������
template<typename T>
vector<vector<T>> mult(const matrix<T>& A1, const matrix<T>& A2) {
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


// ���������������� ���������� ������� (�������� ������� ������� �� ����������������)
template<typename T>
void transpose(vector<vector<T>>& A) {
	size_t n = A.size();
	// ��� ���� ������
	for (size_t i = 0; i < n - 1; ++i)
		for (size_t j = i + 1; j < n; ++j)
			swap(A[i][j], A[j][i]);
}



//����� ������� �
template<typename T>
void print_A(const matrix<T>& vec_A)
{
	for (auto& i : vec_A) {
		for (auto j : i)
			cout << fixed << setprecision(4) << j << "\t";
		cout << endl;
	}
}


//����� �������
template<typename T>
void print_vec(const vector<T>& vec)
{

	//copy(vec.begin(), vec.end(), std::ostream_iterator<double>(std::cout, " "));
	for (size_t i = 0; i < vec.size(); i++)
		cout << fixed << setprecision(8) << vec[i] << "\t";
	cout << endl << endl;
}

//���������� �� �����
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

// ��������� �� ������� �������� T_ij �����
// initIndex ��������� �� ��, � ������ ������� ��������� �������� (����� ��� ������� R � QR-����������)
template<typename T>
void multT(vector<vector<T>>& A, size_t n, size_t ti, size_t tj, T c, T s, size_t initIndex = 0) {
	//size_t n = A.size();
	// ��������, ������ ��� � ����� ��� ������ ����� ��������
	vector<T> Ati(A[ti]);
	vector<T> Atj(A[tj]);

	for (size_t j = initIndex; j < n; ++j) {
		A[ti][j] = c * Ati[j] + s * Atj[j];
		A[tj][j] = -s * Ati[j] + c * Atj[j];
	}
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
			multT(Q, n, i, j, c_ij, s_ij);
			multT(A, n, i, j, c_ij, s_ij, i);
		}
	}
	// ��� �������
	// �������� �� �������������
	for (size_t i = 0; i < n; ++i)
		if (abs(A[i][i]) <= numeric_limits<T>::epsilon() * 10)
			throw "����������� �������";
}

// ����� (���������) �������� ���� ��������
template<typename T>
T norm(const vector<T>& x1, const vector<T>& x2) {
	T res = 0;
	size_t n = x1.size();
	for (size_t i = 0; i < n; ++i)
		res += (x1[i] - x2[i]) * (x1[i] - x2[i]);
	return sqrt(res);
}

// ����� (���������) �������
template<typename T>
T norm(const vector<T>& x) {
	T res = 0;
	size_t n = x.size();
	for (size_t i = 0; i < n; ++i)
		res += x[i] * x[i];
	return sqrt(res);
}

// ��������� ������
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

	// ��������� �����������
	vector<T> xk(n, 0);
	xk[0] = 1;

	//print_vec(xk);

	// ������� ��� ������� ����
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

		// ����� ���������� ����� ��������, ��� �������� ������������ ���� �
		// ��� �� ������, ������ ����������� �� -1
		for (size_t i = 0; i < n; ++i)
			xk[i] *= -1;
		coef = min(coef, norm(xk1, xk));

		xk = xk1;

		//print_vec(xk);
	} while (coef >= eps);

	return xk;
}

// �������� ������� �����������
//���� ������� �����������, ����� ��������������
template<typename T>
void Hessenberg(matrix<T>& A)
{
	// ��� ���� ������� ������� A ������ ���� ����������
	size_t n = A.size();

	for (size_t i = 1; i < n - 1; ++i)
	{
		for (size_t j = i + 1; j < n; ++j) {
			T c = A[i][i - 1];
			T s = A[j][i - 1];

			// ����������
			T denom = sqrt(c * c + s * s);

			//if (denom <= numeric_limits<T>::epsilon() * 10)
			//	continue;
			T alpha = c / denom;
			T beta = s / denom;

			//��������� ������� � ����� �� T_kl
			//for (size_t k = 0; k < n; k++) {
			for (size_t k = i - 1; k < n; k++) {

				// ��������, ������ ��� ��� �������� ����� ��������
				T Aik = A[i][k];
				T Ajk = A[j][k];

				A[i][k] = alpha * Aik + beta * Ajk;
				A[j][k] = -beta * Aik + alpha * Ajk;
			}


			//��������� ������ �� (T_kl)^T ������������ �������
			for (size_t k = 0; k < n; k++) {
				//��������� ��������, ��� ��������
				T Aki(A[k][i]);
				T Akj(A[k][j]);

				A[k][i] = alpha * Aki + beta * Akj;
				A[k][j] = -beta * Aki + alpha * Akj;
			}

		}
	}

}



// ���� RQ ��� ������� �����������. ������ ������� �������
// n - ����� ����������� ��� ������� A ������� ����� � ����
// � �������� ������ ����������� ����� n ����� ����������
template<typename T>
void qrHessenberg(matrix<T>& A, size_t n) {
	for (size_t i = 0; i < n - 1; ++i) {
		size_t j = i + 1;

		T c = A[i][i];
		T s = A[j][i];

		// ����������
		T denom = sqrt(c * c + s * s);
		T alpha = c / denom;
		T beta = s / denom;

		// ��������, ������ ��� � ����� ��� ������ ����� ��������
		vector<T> Aik(A[i]);
		vector<T> Ajk(A[j]);
		//��������� ������� � ����� �� T_kl
		//for (size_t k = i; k < n; k++) {
		// ����� ����� ������, ��� 
		for (size_t k = min(i - 1, i); k < n; k++) {
			A[i][k] = alpha * Aik[k] + beta * Ajk[k];
			A[j][k] = -beta * Aik[k] + alpha * Ajk[k];
		}

		//��������� ������ �� (T_kl)^T ������������ �������
		// TODO �� ����� ��������� ���������� �������� ������� n �� ...
		for (size_t k = 0; k < n; k++) {
			//��������� ��������, ��� ��������
			T aki(A[k][i]);
			T akj(A[k][j]);

			A[k][i] = alpha * A[k][i] + beta * A[k][j];
			A[k][j] = -beta * aki + alpha * akj;
		}
	}
}



// ����� ����������� ����� ������� QR-����������
// ���������� ������ ����������� �����
// TODO ������� ���������� ��������
template<typename T>
vector<T> eigValQR(const matrix<T>& A, T eps) {
	matrix<T> H(A);
	Hessenberg(H);
	size_t n = H.size();
	//cout << "����� ������� ����������� ����� QR-�����������" << endl;
	//print_A(H);
	//cout << endl;

	vector<T> res;
	res.reserve(n);

	size_t iterCount = 0;
	size_t operCount = 0;

	while (n != 1) {
		while (abs(H[n - 1][n - 2]) >= eps) {
			//�������� ��� ������
			T sigma = H[n - 1][n - 1];
			//sigma = 0.0;
			//��������� ����� (-�����)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] -= sigma;
			}
			//������������ ������� �� A_k+1 = Q * A_k * (Q)^T 
			qrHessenberg(H, n);
			operCount += 6 * n * n + 7 * n - 13;
			//��������� �������� ����� (+�����)
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
	res.push_back(H[0][0]); //lambda 1 (��� ������)

	cout << "���-�� ��������: " << iterCount << endl;
	cout << "���-�� ��������: " << operCount << endl;
	return res;
}



// ���� RQ ��� ������� ������� (��/��� ������). ������ ������� �������
// n - ����� ����������� ��� ������� A ������� ����� � ����
// � �������� ������ ����������� ����� n ����� ����������
template<typename T>
void qrHessenbergSimple(matrix<T>& A, size_t n) {
	matrix<T> Q = eig<T>(n);
	for (size_t i = 0; i < n - 1; ++i) {
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
			//�������� ��� ������
			T sigma = H[n - 1][n - 1];
			//sigma = 0.0;
			//��������� ����� (-�����)
			for (size_t i = 0; i < n; ++i) {
				H[i][i] -= sigma;
			}
			//������������ ������� �� A_k+1 = Q * A_k * (Q)^T 
			qrHessenbergSimple(H, n);

			operCount += 10.0 * n * n*n/3.0 + n * n/2.0 - 23*n/6.0;
			//��������� �������� ����� (+�����)
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
	res.push_back(H[0][0]); //lambda 1 (��� ������)

	cout << "���-�� ��������: " << iterCount << endl;
	cout << "���-�� ��������: " << operCount << endl;
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




// ����� ����������� ����� � ����������� �������� � �������������� QR-����������
// � ������ �������� ��������.
// ���������� ������ ��� "����������� �����" - "��������������� ��� �����������
// ������".
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

// ���������� ��������� ����������� � ������������ ������� ���������� �������
template<typename T>
vector<T> generateX0(size_t n) {
	vector<T> res;
	res.reserve(n);
	for (size_t i = 0; i < n; ++i)
		res.push_back(rand());
	normalize(res);
	return res;
}

// ����� ����������� ����� � �������� ��� ������ ����������� �����
// ���������� ������ ��� "������. �����"-"������. ������"
template<typename T>
vector<pair<T, vector<T>>> RayleighRelation(const matrix<T>& A, T eps) {
	size_t n = A.size();
	vector<pair<T, vector<T>>> res;
	res.reserve(n);

	matrix<T> A1(A);

	for (size_t i = 0; i < n; ++i) {
		// ��������� ����������� ������������ �������
		vector<T> xk = generateX0<T>(n);
		// ��� ������� ���. ������� ������������
		// ��� ���������, ����� �� �������� ������ � ���� �� ���. �������
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

			// ������� ��� ������� ����
			matrix<T> M(A1);
			for (size_t i = 0; i < n; ++i)
				M[i][i] -= lamk;

			vector<T> xkCopy(xk);
			xk1 = gauss(M, xkCopy);

			// ������� �� xk1 ���������� ��� ��������� ����������� ��������
			// � ������� ������������� ���� �����
			for (size_t j = 0; j < res.size(); ++j) {
				vector<T> ej = res[j].second;
				T f = mult(xk1, ej);
				for (size_t k = 0; k < n; ++k)
					xk1[k] -= f * ej[k];
			}
			////////////////

			normalize(xk1);

			//�������� �������� ����� �������� ��������� ����������� � ������ ������ �
			//���� �� ������� � ������� �������
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

// ���������, ��� ����� eigVal � ������ eigVec ������������� ��������
// ������������. ���������, ��� (A - eigVal * E) * eigVec = 0.
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
// "������������" ���������� � ������� �����������
// �������� �������� ������� � ������������������ ����
// (�� ������ �� ���������)
// ��� 4 ������������ �������
template<typename T>
void HessenbergIncorr(const matrix<T>& M)
{
	matrix<T> A(M);

	// ��� ���� ������� ������� A ������ ���� ����������
	size_t n = A.size();

	//for (size_t i = 1; i < n - 1; ++i)
	for (size_t i = 0; i < n - 1; ++i)
	{
		for (size_t j = i + 1; j < n; ++j) {
			T c = A[i][i];
			T s = A[j][i];

			// ����������
			T denom = sqrt(c * c + s * s);

			T alpha = c / denom;
			T beta = s / denom;

			//��������� ������� � ����� �� T_kl
			for (size_t k = 0; k < n; k++) {
				// ���� ���������, ������� ������ ������ ��������������
				//for (size_t k = i; k < n; k++) {

					// ��������, ������ ��� ��� �������� ����� ��������
				T Aik = A[i][k];
				T Ajk = A[j][k];

				A[i][k] = alpha * Aik + beta * Ajk;
				A[j][k] = -beta * Aik + alpha * Ajk;
			}


			//��������� ������ �� (T_kl)^T ������������ �������
			for (size_t k = 0; k < n; k++) {
				//��������� ��������, ��� ��������
				T Aki(A[k][i]);
				T Akj(A[k][j]);

				A[k][i] = alpha * Aki + beta * Akj;
				A[k][j] = -beta * Aki + alpha * Akj;
			}
			//print_A(A);

		}
	}
	cout << "\"������������\" ����������:" << endl;
	print_A(A);
	cout << "\"����������\" ����������:" << endl;
	matrix<T> A2(M);
	Hessenberg(A2);
	print_A(A2);
}



// ����� �������� �������� � ������� ���������� �����������
template<typename T>
vector<T> eigenvector(const matrix<T> A, T eigVal, const vector<T>& x0, T eps) {
	size_t n = A.size();

	// ��������� �����������
	vector<T> xk = x0;

	// ������� ��� ������� ����
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

		// ����� ���������� ����� ��������, ��� �������� ������������ ���� �
		// ��� �� ������, ������ ����������� �� -1
		for (size_t i = 0; i < n; ++i)
			xk[i] *= -1;
		coef = min(coef, norm(xk1, xk));

		xk = xk1;
	} while (coef >= eps);

	return xk;
}


// ����� ������ ������������ ����� � ������� ��� ������ ����������� �����
// ���������� ���� "������. �����"-"������. ������"
// ��������� �� ���� ��������� ����������� x0
template<typename T>
pair<T, vector<T>> RayleighRelationSimple(const matrix<T>& A, const vector<T>& x0, T eps) {
	size_t n = A.size();

	matrix<T> A1(A);


	// ��������� �����������
	vector<T> xk = x0;
	normalize(xk);

	T lamk;

	vector<T> xk1;

	T coef;
	do {
		lamk = mult(mult(A1, xk), xk);

		// ������� ��� ������� ����
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

// ������ � 6 ������������ �������
void question6() {
	matrix<double> A;
	size_t N = out_of_file(A, "EIGEN1.TXT");

	double eps = 1e-3;

	double lam1 = 57.165516324888;
	// ����������� ������ � ����������� ������ 57.165516324888
	vector<double> x1{ -0.063563494887, 0.443928222700, 0.888001562853, -0.101688935383 };

	double lam2 = 22.359615081837;
	// ����������� ������ � ����������� ������ 22.359615081837
	vector<double> x2{ 0.777968847720, 0.582505207755, -0.229455791959, 0.052935757783 };

	cout << "����� �������� ��������:";
	// �������� � x1 (��-�� lam1)
	auto eInv = eigenvector(A, lam1, x2, eps);
	print_vec(eInv);
	cout << "����������� �����:";
	// �������� � x2 (��-�� ���������� �����������)
	auto lamERayleigh = RayleighRelationSimple(A, x2, eps);
	auto eRayleigh = lamERayleigh.second;
	print_vec(eRayleigh);
}



//����� QR-����������
//R = P^(-1) * A * P
//P - �������������
//
//
//QR - ���������� ������� ������� (O(n3)) ��������
//������ ����� �������� QR-���������� -> �������� ������� � � ����� �����������

int main()
{
	setlocale(LC_ALL, "Russian");
	cout << fixed;
	cout.precision(10);


	typedef double MyType;
	MyType eps = 1e-3;
	matrix<MyType> A;
	matrix<MyType> Q;

	//������ �� �����
	//size_t N = out_of_file(A, "���� 1.txt");
	//size_t N = out_of_file(A, "���� 2.txt"); //�������������� �������
	//size_t N = out_of_file(A, "EIGEN1.TXT");
	size_t N = out_of_file(A, "EIGEN20.TXT");



	//matrix<MyType> A_H(A);
	//cout << "������� �:" << endl;
	//print_A(A);
	//cout << endl << "����� ������� �����������" << endl;
	//Hessenberg(A_H);
	//cout << endl;
	//qrHessenberg(A_H, A_H.size());
	//print_A(A_H);

	//transpose(A);
	////QR-����������
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\tQR - ����������" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	auto eigsQR = eigValVecQR(A, eps);
	for (auto& eig : eigsQR) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "����������� ����� " << fixed << setprecision(10) << eigVal << endl << "����������� ������:";
		print_vec(eigVec);
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}


	/*QR-���������� (��� ������� �����������)*/
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\tQR - ���������� (��� ������� �����������)" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;

	auto eigsQRSimple = eigValVecQRSimple(A, eps);
	for (auto& eig : eigsQRSimple) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "����������� ����� " << fixed << setprecision(10) << eigVal << endl << "����������� ������:";
		print_vec(eigVec);
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}




	//����������� �����
	cout << endl << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	cout << "\t\t\t����������� �����" << endl;
	cout << "----------------------------------------------------------------------------------" << endl;
	//srand(time(NULL));
	auto eigs = RayleighRelation(A, eps);
	for (auto& eig : eigs) {
		MyType eigVal = eig.first;
		vector<MyType> eigVec = eig.second;

		cout << "����������� ����� " << fixed << setprecision(10) << eigVal << endl << "����������� ������:";
		print_vec(eigVec);
		cout << endl;
		bool check = checkEigValAndVec(A, eigVal, eigVec, eps);
	}


	question6();
}



//����� �� 4 ������, ��� � �����, ������������.����� ����� ������ ����,
//��� ����� �������� �� �� �������� � ������������ ����.�� ������� ����, 
//����� � �������� ��� ����, ����� ������, ��� ���� �������.
//��� 6 ������� ������� ������.
//��� � ������� � ����������� �������� ���������� �������� � �������� ���
//QR - ��������� ��� ���������� ����� �����������. � ����� ��������, ��� 
//��������� ����������� ����� � ������� ��� ����������������� �������
//(��� �������������� �������, �������).�� � �������� ������
//��� ����� ����� ��������� ��� ����������� ������� ����������������� �������, 
//�� ����� - �� ����� ���������� ��� ����������� ������������ � ������� ����� 
//� ������������� �����������.�� �����, � �� �����