#pragma once
#include <iostream>
#include <fstream>
#include <complex>
using namespace std;

template<typename _Type>
void CreateMatrix(int _I, int _J, _Type**& _A)
{
    int i1, i2;
    _A = new _Type *[_I];
    for (i1 = 0; i1 < _I; i1++) {
        _A[i1] = new _Type[_J];
        for (i2 = 0; i2 < _J; i2++) {
            _A[i1][i2] = 0.0;
        }
    }
}

template<typename _Type>
_Type** CreateMatrix(int _I, int _J)
{
    _Type** A;
    CreateMatrix(_I, _J, A);
    return A;
}

template<typename _Type>
void CreateMatrix(int _I, _Type**& _A)
{
    CreateMatrix(_I, _I, _A);
}

template<typename _Type>
void CreateVec(int _I, _Type*& vec) {
    vec = new _Type[_I];
}

template<typename _Type>
_Type* CreateVec(int _I) {
    _Type* vec;
    CreateVec(_I, vec);
    return vec;
}

template<typename _Type>
void deleteMemory(int _I, _Type **& _A) {
    int i1;
    for (i1 = 0; i1 < _I; i1++) {
        delete _A[i1];
    }
    delete[]_A;
}

template<typename _Type>
void deleteMemory(_Type*& _A) {
    delete[]_A;
}

template<typename _Type>
double normaVec(_Type* _vec_1, complex<double>* _vec_2, int _N) {
    double norma = 0;
    for (size_t i = 0; i < _N; i++)
    {
        norma += pow(abs(_vec_1[i] - _vec_2[i]), 2);
    }
    return sqrt(norma);
}

template<typename _Type>
double normaVec(_Type* _vec, int _N) {
    double norma = 0;
    for (size_t i = 0; i < _N; i++)
    {
        norma += pow(abs(_vec[i]), 2);
    }
    return sqrt(norma);
}

template<typename _Type>
_Type scalar_mult(_Type* _vec_1, _Type* _vec_2, int _N) {
    _Type scalar = 0;
    for (size_t i = 0; i < _N; i++)
    {
        scalar += _vec_1[i] * _vec_2[i];
    }
    return scalar;
}

template<typename _Type>
_Type scalar_mult_conj(_Type* _vec_conj_1, _Type* _vec_2, int _N) {
    _Type scalar = 0;
    for (size_t i = 0; i < _N; i++)
    {
        scalar += conj(_vec_conj_1[i]) * _vec_2[i];
    }
    return scalar;
}

template<typename _Type>
void copy_vec(_Type* _from_vec, _Type*& _to_vec, int _N) {
    for (size_t i = 0; i < _N; i++)
    {
        _to_vec[i] = _from_vec[i];
    }
}

template<typename _Type>
_Type* copy_vec(_Type* _from_vec, int _N) {
    _Type* new_vec = CreateVec<_Type>(_N);
    for (size_t i = 0; i < _N; i++)
    {
        new_vec[i] = _from_vec[i];
    }
    return new_vec;
}

template<typename _Type>
void minus_vec(_Type* _vec_1, _Type* _vec_2, int _N, _Type*& _result, _Type k_1 = 1, _Type k_2 = 1) {
    for (size_t i = 0; i < _N; i++)
    {
        _result[i] = _vec_1[i] * k_1 - _vec_2[i] * k_2;
    }
}

template<typename _Type>
void plus_vec(_Type* _vec_1, _Type* _vec_2, int _N, _Type*& _result, _Type k_1 = 1, _Type k_2 = 1) {
    for (size_t i = 0; i < _N; i++)
    {
        _result[i] = _vec_1[i] * k_1 + _vec_2[i] * k_2;
    }
}


void mult_matr_vec(complex<double>** _matr, complex<double>* _vec, int _N, complex<double>*& result) {

    for (size_t i = 0; i < _N; i++)
    {
        result[i] = 0;
        for (size_t j = 0; j < _N; j++)
        {
            result[i] += _matr[i][j] * _vec[j];
        }
    }
}

complex<double>* mult_matr_vec(complex<double>** _matr, complex<double>* _vec, int _N) {
    complex<double>* result = new complex<double>[_N];

    for (size_t i = 0; i < _N; i++)
    {
        result[i] = 0;
        for (size_t j = 0; j < _N; j++)
        {
            result[i] += _matr[i][j] * _vec[j];
        }
    }
    return result;
}


void printMatrix(complex <double>** A, int N, int M) {
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < M; j++) {
            cout << A[i][j] << "  ";
        }
        cout << endl;
    }
}

void mult_matr_matr(complex<double>** _matr_1, complex<double>** _matr_2, int _N, complex<double>**& result) {
#pragma omp parallel
    {
#pragma omp for 
        for (int i = 0; i < _N; i++)
        {

            for (size_t j = 0; j < _N; j++)
            {
                result[i][j] = 0;
                for (size_t k = 0; k < _N; k++)
                {
                    result[i][j] += _matr_1[i][k] * _matr_2[k][j];
                }
                //result[i] += _matr[i][j] * _vec[j];
            }
            //cout << "mult I = " << i << endl;
        }
    }
}

complex<double>** mult_matr_matr(complex<double>** _matr_1, complex<double>** _matr_2, int _N) {
    complex<double>** result;
    CreateMatrix(_N, _N, result);
    mult_matr_matr(_matr_1, _matr_2, _N, result);
    return result;
}

template<typename type_ptr>
void switchPTR(type_ptr*& first_ptr, type_ptr*& second_ptr) {
    type_ptr* tmp = first_ptr;
    first_ptr = second_ptr;
    second_ptr = tmp;
}

complex<double>* BICGstab(complex<double>** _matr, complex<double>* _vec, int _N) {
    complex<double>* r_k = new complex<double>[_N];
    complex<double>* r_volna = new complex<double>[_N];
    complex<double>* x_km1 = new complex<double>[_N];
    complex<double>* x_k = new complex<double>[_N];

    complex<double>* tmp = new complex<double>[_N];

    complex<double> rho, rho_km1, alpha, beta, omega;

    complex<double>* v_k = new complex<double>[_N];
    complex<double>* p_k = new complex<double>[_N];
    complex<double>* s = new complex<double>[_N];
    complex<double>* t = new complex<double>[_N];

    copy_vec(_vec, x_km1, _N);

    mult_matr_vec(_matr, x_km1, _N, tmp);
    minus_vec(_vec, tmp, _N, r_k);
    copy_vec(r_k, r_volna, _N);

    rho = 1;
    rho_km1 = 1;
    alpha = 1;
    beta = 1;
    omega = 1;

    for (size_t i = 0; i < _N; i++)
    {
        v_k[i] = 0;
        p_k[i] = 0;
    }
    int iter = 0;
    do {
        /////1
        rho = scalar_mult_conj(r_volna, r_k, _N);
        /////1
        /////2
        beta = rho / rho_km1 * alpha / omega;
        /////2
        /////3
        minus_vec(p_k, v_k, _N, tmp, complex<double>(1, 0), omega);
        plus_vec(r_k, tmp, _N, p_k, complex<double>(1, 0), beta);
        /////3
        /////4
        mult_matr_vec(_matr, p_k, _N, v_k);
        /////4
        /////5
        alpha = rho / scalar_mult_conj(r_volna, v_k, _N);
        /////5
        /////6
        minus_vec(r_k, v_k, _N, s, complex<double>(1, 0), alpha);
        /////6
        /////7
        mult_matr_vec(_matr, s, _N, t);
        /////7
        /////8
        omega = scalar_mult(t, s, _N) / scalar_mult(t, t, _N);
        /////8
        /////9
        plus_vec(s, p_k, _N, tmp, omega, alpha);
        plus_vec(x_km1, tmp, _N, x_k);
        /////9
        /////10
        minus_vec(s, t, _N, r_k, complex<double>(1, 0), omega);
        /////10

        rho_km1 = rho;
        switchPTR(x_k, x_km1);
        cout << "norma b = " << normaVec(x_k, x_km1, _N) << endl;
        cout << "norma r = " << normaVec(r_k, _N) << endl;
        //cout << "norma b = " << normaVec(x_k, x_km1, _N) << endl;
        iter++;
    } while (normaVec(x_k, x_km1, _N) > 0.00000001 && normaVec(r_k, _N) > 0.000001);
    cout << "Iter = " << iter << endl;
    return x_km1;
}

void reverse_triangle_down(complex<double>** _matr, int _N, complex<double>** result) {
    /*complex<double>** obr;
    CreateMatrix(_N, _N, obr);*/

    /*for (size_t i = 0; i < _N; i++)
    {
        
        for (size_t j = 0; j < _N; j++)
        {
            if (i == j)
                result[i][j] = 1. / _matr[i][j];
            else if (j < i) {
                result[i][j] = 0;
                for (size_t k = 0; k <= i; k++)
                {
                    result[i][j] += _matr[i][k] * _matr[j][k];
                }
                result[i][j] /= _matr[j][j];
            }
            else
                result[i][j] = 0.;
        }
    }*/
    for (size_t i = 0; i < _N; i++)
    {
        result[i][i] = 1. / _matr[i][i];
    }
    for (size_t j = 0; j < _N - 1; j++)
    {
        for (size_t i = j + 1; i < _N; i++)
        {
            //cout << "i = " << i << endl;
            result[i][j] = 0;
            for (size_t k = j; k < i; k++)
            {
                //cout << _matr[i][k] << " * " << result[k][j] << endl;
                result[i][j] += _matr[i][k] * result[k][j];
            }
            result[i][j] *= (-1. / _matr[i][i]);
        }
    }

}

void Gauss(complex <double>** Matrix, complex <double>* Vec, int Nm) {
    complex <double> ed(1.0, 0.0);
    complex <double> nul(0.0, 0.0);

    for (int k = 0; k < Nm; k++) {

        cout << "k = " << k << "    " << endl;

        if (Matrix[k][k] != ed) {
            complex <double> T = Matrix[k][k];
            for (int j = k; j < Nm; j++) {//нормирование строки
                Matrix[k][j] = Matrix[k][j] / T;
            }
            Vec[k] = Vec[k] / T;
        }
        for (int i = k; i < Nm; i++) { //проходим по столбцу
            if ((Matrix[i][k] != ed) && (i != k)) {
                complex <double> T = Matrix[i][k];
                Matrix[i][k] = 0;
                for (int j = k + 1; j < Nm; j++) { //проходим по двум строкам и вычитаем их
                    Matrix[i][j] -= Matrix[k][j] * T;
                }
                Vec[i] -= Vec[k] * T;
            }
        }

    }
    for (int i = Nm - 1; i >= 0; i--) {
        complex <double> Sum = Vec[i];
        for (size_t j = i + 1; j < Nm; j++) {
            Sum -= Matrix[i][j] * Vec[j];
        }
        Vec[i] = Sum;
    }
}

template<typename type_ptr>
double norma_matr(type_ptr** _matr, int _N) {
    double norma = 0;
    for (size_t i = 0; i < _N; i++)
        for (size_t j = 0; j < _N; j++)
            norma += abs(_matr[i][j]);
    return norma;
}

void alghoritm_hottenlinga(complex <double>** _Matrix, complex <double>** _begin_m, int _N, complex <double>**result) {
    complex <double>** tmp;
    complex <double>** D_kp1;
    complex <double>** R;

    CreateMatrix(_N, _N, tmp);
    CreateMatrix(_N, _N, D_kp1);
    CreateMatrix(_N, _N, R);

    
    for (size_t iter = 0; iter < 40; iter++)
    {
        mult_matr_matr(_Matrix, _begin_m, _N, R);
        for (size_t i = 0; i < _N; i++)
            for (size_t j = 0; j < _N; j++)
                R[i][j] = (i == j) ? (1.0 - R[i][j]) : (-R[i][j]);

        cout << "norma alg = " << norma_matr(R, _N) << endl;

        for (size_t i = 0; i < _N; i++)
            for (size_t j = 0; j < _N; j++)
                R[i][j] = (i == j) ? (1.0 + R[i][j]) : (R[i][j]);

        mult_matr_matr(_begin_m, R, _N, D_kp1);
        switchPTR(_begin_m, D_kp1);
    }
    


}


complex<double>** predobusl(complex<double>** _Matrix, int N, double omega = 1.0) {
    complex<double>** D_0;
    complex<double>** D_12;
    complex<double>** tmp;
    complex<double>** L;
    CreateMatrix(N, N, D_0);
    CreateMatrix(N, N, D_12);
    CreateMatrix(N, N, tmp);
    CreateMatrix(N, N, L);

    for (size_t i = 1; i < N; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            L[i][j] = _Matrix[i][j];
        }
    }
    //printMatrix(_Matrix, N, N);
    //printMatrix(L, N, N);
    //cout << endl;
    for (size_t i = 0; i < N; i++)
    {
        D_12[i][i] = sqrt(2 - omega) / sqrt((1. / omega) * _Matrix[i][i]);
    }
    for (size_t i = 0; i < N; i++)
    {
        D_0[i][i] = 1.0 / ((1. / omega) * _Matrix[i][i]);
        //cout << " i = " << i << "   el = " << D_0[i][i] << endl;
    }

    mult_matr_matr(L, D_0, N, tmp);

    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            tmp[i][j] = (i == j) ? (1.0 - tmp[i][j]) : (-tmp[i][j]);
        }

    }

    mult_matr_matr(D_12, tmp, N, D_0);
    //return D_0;

    //printMatr(D_0, N);
    //exit(-1);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < i; j++)
        {
            tmp[i][j] = D_0[j][i];
        }
    }

    mult_matr_matr(tmp, D_0, N, D_12);      //транспонированная на обычную
    /*
    N = 4;
    D_12[0][0] = 3;
    D_12[1][0] = 1;
    D_12[2][0] = 2;
    D_12[3][0] = 1;
    D_12[1][1] = 8;
    D_12[2][1] = 3;
    D_12[3][1] = 2;
    D_12[2][2] = 6;
    D_12[3][2] = 4;
    D_12[3][3] = 4;*/
    //printMatr(D_12, N);

    //reverse_triangle_down(D_12, N, tmp);
    //printMatr(tmp, N);
    //return tmp;

    //mult_matr_matr(tmp, D_12, N, D_0);

    //cout << "___________\n";
    //printMatr(D_0, N);
    //cout << "___________\n";

    //exit(-1);

    /*for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            D_0[i][j] = E[i][j] - tmp[i][j];
        }
    }
    double norma = 0;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            norma += abs(D_0[i][j]);
        }

    }*/
    //cout << "abs D_0 = " << sqrt(norma) << endl;

    //mult_matr_matr(_Matrix, D_12, N, tmp);

    //for (size_t i = 0; i < N; i++)
    //{
    //    for (size_t j = 0; j < N; j++)
    //    {
    //        if (i == j) {
    //            tmp[i][j] = 1. - tmp[i][j];
    //        }
    //        else {
    //            tmp[i][j] = -tmp[i][j];
    //        }
    //    }

    //}
    //double norma = 0;
    //for (size_t i = 0; i < N; i++)
    //{
    //    for (size_t j = 0; j < N; j++)
    //    {
    //        norma += abs(tmp[i][j]);
    //    }
    //}
    //cout << "norma = " << norma << endl;

    alghoritm_hottenlinga(_Matrix, D_12, N, tmp);
    //exit(-1);

    return D_12;

}


void MatrixOnVector(int n, int flag, complex<double>** A, complex<double>* B, complex<double>* W, complex<double>* U)
{
    int i, j;

    // умножение обычной матрицы на вектор
    if (1 == flag) {
        for (i = 0; i < n; i++) {
            U[i] = 0.0;
            for (j = 0; j < n; j++) {
                U[i] += A[i][j] * B[j]; // правая часть
            }
        }
    }
    else {// умножение сопряженной матрицы на вектор
        for (i = 0; i < n; i++) {
            U[i] = 0.0;
            for (j = 0; j < n; j++) {
                U[i] += conj(A[j][i]) * B[j]; // правая часть
            }
        }
    }

    for (i = 0; i < n; i++) {
        U[i] *= W[i];
    }

}
complex<double> Norma2(int n, complex<double>* x)
{
    int i;
    complex<double> s, ret;

    complex<double>* xs;
    xs = new complex<double>[n];

    for (i = 0; i < n; i++) {
        xs[i] = real(x[i]) - complex<double>(0, 1.) * imag(x[i]);
    }

    s = 0.0;
    for (i = 0; i < n; i++) {
        s = s + x[i] * xs[i];
    }

    ret = s;

    delete[]xs;

    return ret;
}

int Gradient(int n, complex<double>** A, complex<double>* V, complex<double>* W, complex<double>* U)
{
    int i;				//вспомагательные переменые

    int k;				// число итераций
    complex<double> b, s, al;		//вспомагательные переменые
    complex<double>	h2, Nf02;	//вспомагательные переменые


    //выделение памяти под вектора
    complex<double>* f = CreateVec<complex<double>>(n);
    complex<double>* f0 = CreateVec<complex<double>>(n);
    complex<double>* w = CreateVec<complex<double>>(n);
    complex<double>* q = CreateVec<complex<double>>(n);
    complex<double>* h = CreateVec<complex<double>>(n);
    complex<double>* p = CreateVec<complex<double>>(n);


    // умножаем на вектор геометрии
    for (i = 0; i < n; i++) {
        V[i] *= W[i];
    }

    //начальное приближение
    for (i = 0; i < n; i++) {
        f[i] = V[i];
        U[i] = V[i];
    }

    // правая часть нового матричного уравнения А1*Аu = А1*f  и векторы w = A*u и q = A1*w
    MatrixOnVector(n, 0, A, f, W, f0);
    MatrixOnVector(n, 1, A, U, W, w);
    MatrixOnVector(n, 0, A, w, W, q);

    for (i = 0; i < n; i++) {
        h[i] = q[i] - f0[i];
        p[i] = h[i];
    }

    Nf02 = Norma2(n, f0);	// скалярное произведение
    h2 = Norma2(n, h);	// скалярное произведение

    k = 0;
    const double eps = 0.00000000000000001;	// квадрат погрешности
    // запуск итераций
    while (eps < abs(h2 / Nf02))
    {
        cout << "Norma = " << abs(h2 / Nf02) << endl;
        //cout<<abs(h2/Nf02)<<endl;
        MatrixOnVector(n, 1, A, p, W, w);
        MatrixOnVector(n, 0, A, w, W, q);

        al = h2 / Norma2(n, w);


        //новое приближение решения
        for (i = 0; i < n; i++) {
            U[i] = U[i] - al * p[i];
            h[i] = h[i] - al * q[i];
        }

        b = 1.0 / h2;
        h2 = Norma2(n, h);
        b = b * h2;

        //новое значение p
        for (i = 0; i < n; i++) {
            p[i] = h[i] + b * p[i];
        }

        k++;
        cout << "Num iter = " << k << endl;
    }

    



    return 1;
}