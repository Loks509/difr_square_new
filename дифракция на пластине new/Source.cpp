#include <iostream>
#include <fstream>
#include <complex>
#include "matrix_func.h"
#include <mpi.h>

using namespace std;

const double Pi = acos(-1);

const double v_c = 299792458;
const double Kilo = pow(10, 3);
const double Mega = pow(10, 9);
const double Tera = pow(10, 12);

double freq = 10.0 * Mega;


int start_index = 37;    //менять это!
const int c_of_iter = 100;  //это по желанию


//с этми значениями решается 24 на 24 с трудом!
//double freq = 10008.0 * Mega;
//int N_int = 30;
//const int n = 24;

double K0 = 2. * Pi * freq / v_c;

int N_int = 24;

const int n = 16;
const int N = n * n;
const double  A = -1., B = 1.,
C = -1., D = 1.;


const int n_obr = n;
const int N_obr = n_obr * n_obr;

const double min_good = 50;
const double max_good = 150;

const double min_bad = 150;
const double max_bad = 1581;

double x_ist = -2, y_ist = 2;

int count_mass[n * n] = {0};

int global_rank = -1;

complex<double> _H(double x) {
    return complex<double>(_j0(x), _y0(x));
}

int get_random_index() {
    return rand() % 16;
}

double K_f_1(double x, double y) {
    static int* var = nullptr;
    if (var == nullptr) {
        var = new int[20];
        cout << get_random_index() << endl;
        var[0] = get_random_index();
        var[1] = get_random_index();
        var[2] = get_random_index();
        var[3] = get_random_index();
        var[4] = get_random_index();
        var[5] = get_random_index();
        var[6] = get_random_index();
        var[7] = get_random_index();
        var[8] = get_random_index();
        var[9] = get_random_index();
        var[10] = get_random_index();
        var[11] = get_random_index();
        var[12] = get_random_index();
        var[13] = get_random_index();
        var[14] = get_random_index();
        var[15] = get_random_index();
    }

    if (x == var[0] && y == var[1])
        return 300;
    else if ((x == var[2] || x == (var[2] + 1)) && (y == var[3] || y == (var[3] + 1)))
        return 500;
    else if (x == var[4] && y == var[5])
        return 500;
    else if (x == var[6] && y == var[7])
        return 500;
    else if (x == var[8] && y == var[9])
        return 500;
    else if (x == var[10] && y == var[11])
        return 500;
    else if (x == var[12] && y == var[13])
        return 500;
    else if ((x == var[14] || x == (var[14] - 1) || x == (var[14] + 1)) && (y == (var[15]) || y == (var[15] - 1) || y == (var[15] + 1)))
        return 500;
    else
        return 100;
}


double K_f(double x, double y) {
    return K_f_1(x, y);
}
complex<double> Kernel(double k, double x_1, double x_2, double y_1, double y_2) {
    complex <double> ed(0, 1.0);
    double l;
    l = sqrt(pow(x_1 - y_1, 2) + pow(x_2 - y_2, 2));
    return _H(k * l);
    //return exp(ed * K0 * l) / (4.0 * Pi * l);
}

complex<double> Integr(double x_beg, double x_end, double y_beg, double y_end, double x_koll, double y_koll, double K, int I, int J) {
    double h_x = (x_end - x_beg) / (double)N_int;
    double h_y = (y_end - y_beg) / (double)N_int;

    double x_c = (x_end + x_beg) / 2.0; //альтренативный вариант неоднородности
    double y_c = (y_end + y_beg) / 2.0;

    complex<double> Sum = 0;
    for (size_t i = 0; i < N_int; i++)
    {
        double x = x_beg + i * h_x + h_x / 2.;
        for (size_t j = 0; j < N_int; j++)
        {
            double y = y_beg + j * h_y + h_y / 2.;
            Sum += (pow(K_f(I, J), 2) - pow(K0, 2)) * Kernel(K, x, y, x_koll, y_koll);
        }
    }
    return Sum * h_x * h_y;
}


complex<double> Integr_Revers(double x_beg, double x_end, double y_beg, double y_end, double x_koll, double y_koll, double K) {
    double h_x = (x_end - x_beg) / (double)N_int;
    double h_y = (y_end - y_beg) / (double)N_int;
    complex<double> Sum = 0;
    for (size_t i = 0; i < N_int; i++)
    {
        double x = x_beg + i * h_x + h_x / 2.;
        for (size_t j = 0; j < N_int; j++)
        {
            double y = y_beg + j * h_y + h_y / 2.;
            Sum += Kernel(K0, x, y, x_koll, y_koll);
        }
    }
    return Sum * h_x * h_y;
}


//плоская волна
//complex<double> fallWave(double k, double x) {
//    complex <double> ed(0, 1.0);
//    return exp(ed * k * x);
//}

//сферическая волна
complex<double> fallWave(double k, double x, double y, double x_c = -2, double y_c = -2) {
    complex <double> ed(0, 1.0);
    //return exp(ed * k * x);
    double r = sqrt(pow(x - x_c, 2) + pow(y - y_c, 2));
    return ed / 4. * _H(k * r);
    //return exp(ed * k * r);
}

void printInFile(complex<double>* _Vec, double _A, double _C, double _h_x, double _h_y, int _n, string name_file = "default_name.txt") {
    int N = _n * _n;
    ofstream file(name_file);
    file << "X Y Z F R I" << endl;
    for (size_t i = 0; i < N; i++)
    {
        int i_koll = i / _n;
        int j_koll = i % _n;
        double x = _A + i_koll * _h_x + _h_x / 2.;
        double y = _C + j_koll * _h_y + _h_y / 2.;
        file << x << " " << y << " 0 " << abs(_Vec[i]) << " " << _Vec[i].real() << " " << _Vec[i].imag() << endl;
    }
    file.close();
}

complex<double> getIntensivity(double _x, double _y, double _A, double _C, double _h_x, double _h_y, int _n, complex<double>* _vec) {
    complex<double> Intens = fallWave(K0, _x, _y, x_ist, y_ist);
    for (size_t i = 0; i < _n * _n; i++)
    {
        int i_int = i / _n;
        int j_int = i % _n;
        double x_beg = _A + i_int * _h_x;
        double y_beg = _C + j_int * _h_y;
        double x_end = x_beg + _h_x;
        double y_end = y_beg + _h_y;
        Intens += Integr(x_beg, x_end, y_beg, y_end, _x, _y, K0, i_int, j_int) * _vec[i];
    }
    return Intens;
}

void print_point_view(double** _points, int _N) {
    ofstream file("points_of_view.txt");
    file << "X Y Z" << endl;
    for (size_t i = 0; i < _N; i++)
    {
        file << _points[i][0] << " " << _points[i][1] << " 0 \n";
    }
}

complex<double> create_noise(complex<double> value, double percent) {
    complex<double> max = value / 100.0 * percent;

    double noise = (rand() % 2000) / 1000.0 - 1.0;
    //cout << noise*max <<"  "<<value<< endl;
    return value + noise * max;
}

complex<double>* filter(complex<double>* Matr, int N, int iter) {
    static complex<double>* Result = nullptr;
    if (Result == nullptr)
        CreateVec(N, Result);
    for (size_t i = 0; i < N; i++)
    {
        if (pow(min_good, 2) < Matr[i].real() && Matr[i].real() < pow(max_good, 2)) {
            Matr[i] = 0;  //фон
            count_mass[i]++;
            Result[i] = abs(Result[i].real()) * (double)(count_mass[i] - 1) + abs(Matr[i].real());
            Result[i] /= count_mass[i];
        }
        else if (pow(min_bad, 2) < Matr[i].real() && Matr[i].real() < pow(max_bad, 2)) {
            count_mass[i]++;
            Result[i] = abs(Result[i].real()) * (double)(count_mass[i] - 1) + abs(Matr[i].real());
            Result[i] /= count_mass[i];
        }
        else
            Matr[i] = 0;
    }
    return Result;
}

int main() {
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    ofstream main_file(to_string(start_index+global_rank)+"_experiment.txt");
    ofstream ans(to_string(start_index + global_rank) + "_answer.txt");
    srand(time(0) + global_rank);
    complex<double>** Am;
    complex<double>* Vec;
    CreateMatrix(N, N, Am);
    CreateVec(N, Vec);

    
    double h_x = (B - A) / (double)n;
    double h_y = (D - C) / (double)n;
    //cout << h_x << " " << h_y << endl;
    for (size_t iteration = 0; iteration < c_of_iter; iteration++)
    {
        if (global_rank == 0) {
            cout << "Current iteration: " << iteration << endl;
        }
        /*if (iteration < c_of_iter) {
            freq += 0.5 * Mega;
        }
        else {
            x_ist =  sqrt(8) * cos(iteration / (double)c_of_iter * 2 * Pi);
            y_ist =  sqrt(8)* sin(iteration / (double)c_of_iter * 2 * Pi);
            freq = 5.5*Mega;
        }*/
        freq += 0.1 * Mega;
        /*x_ist = sqrt(8) * cos(iteration / (double)c_of_iter * 2 * Pi);
        y_ist = sqrt(8) * sin(iteration / (double)c_of_iter * 2 * Pi);*/
        //cout << x_ist << "   " << y_ist << endl;
        K0 = 2. * Pi * freq / v_c;

        for (size_t I = 0; I < N; I++)
        {
            int i_koll = I / n;
            int j_koll = I % n;
            double x_koll = A + i_koll * h_x + h_x / 2.;
            double y_koll = C + j_koll * h_y + h_y / 2.;

            for (size_t J = 0; J < N; J++)
            {
                int i_int = J / n;
                int j_int = J % n;
                double x_beg = A + i_int * h_x;
                double y_beg = C + j_int * h_y;
                double x_end = x_beg + h_x;
                double y_end = y_beg + h_y;
                //cout << x_beg << " " << x_end << " " << y_beg << " " << y_end << endl;
                Am[I][J] = (I == J) ? 1 : 0;
                Am[I][J] -= Integr(x_beg, x_end, y_beg, y_end, x_koll, y_koll, K0, i_int, j_int);
            }
            Vec[I] = fallWave(K0, x_koll, y_koll, x_ist, y_ist);
            //cout << "\rI = " << I;
        }
        //cout << endl;
        Gauss(Am, Vec, N);
        printInFile(Vec, A, C, h_x, h_y, n);
        /*ofstream file("pole.txt");
        file << "X Y Z F R I" << endl;
        for (double x = -2; x <= 2; x += 0.1) {
            cout << "x = " << x << endl;
            for (double y = -2; y <= 2; y += 0.1) {
                auto var = getIntensivity(x, y, A, C, h_x, h_y, n, Vec);
                file << x << " " << y << " 0 " << abs(var) << " " << var.real() << " " << var.imag() << endl;
            }
        }
        file.close();*/

        ///////обратная задача
        //создание точек наблюдения
        double h_x_obr = (B - A) / n_obr;
        double h_y_obr = (D - C) / n_obr;

        double begin_y_pv = 0.1;
        double h_y_pv = .2;
        double** point_view = CreateMatrix<double>(N_obr, 2);

        if (n_obr % 2 != 0) {
            cout << "Error in N_obr\n";
            return -1;
        }

        //for (size_t i = 0; i < n_obr; i++)  //изменение x
        //{
        //    for (size_t j = 0; j < n_obr / 2; j++)    //изменение y
        //    {
        //        point_view[i * n_obr + j][0] = A + i * h_x_obr + h_x_obr / 2.;
        //        point_view[i * n_obr + j][1] = C - begin_y_pv - j * h_y_pv;
        //        point_view[i * n_obr + j + n_obr / 2][0] = A + i * h_x_obr + h_x_obr / 2.;
        //        point_view[i * n_obr + j + n_obr / 2][1] = D + begin_y_pv + j * h_y_pv;
        //    }
        //}

        for (size_t i = 0; i < n_obr; i++)  //изменение x
        {
            for (size_t j = 0; j < n_obr / 4; j++)    //изменение y
            {
                point_view[i * n_obr + j][0] = A + i * h_x_obr + h_x_obr / 2.;
                point_view[i * n_obr + j][1] = C - begin_y_pv - j * h_y_pv;
                point_view[i * n_obr + j + n_obr / 4][0] = A + i * h_x_obr + h_x_obr / 2.;
                point_view[i * n_obr + j + n_obr / 4][1] = D + begin_y_pv + j * h_y_pv;
                point_view[i * n_obr + j + 2 * n_obr / 4][0] = B + j * 0.2;
                point_view[i * n_obr + j + 2 * n_obr / 4][1] = C + i * h_y_obr + h_y_obr / 2.;
                point_view[i * n_obr + j + 3 * n_obr / 4][0] = A - j * 0.2;
                point_view[i * n_obr + j + 3 * n_obr / 4][1] = C + i * h_y_obr + h_y_obr / 2.;
            }
        }
        print_point_view(point_view, N_obr);
        complex<double>** Am_obr, * Vec_obr;
        CreateMatrix(N_obr, Am_obr);
        CreateVec(N_obr, Vec_obr);

        for (size_t I = 0; I < N_obr; I++)
        {
            for (size_t J = 0; J < N_obr; J++)
            {
                int i_int = J / n_obr;
                int j_int = J % n_obr;
                double x_beg = A + i_int * h_x_obr;
                double y_beg = C + j_int * h_y_obr;
                double x_end = x_beg + h_x_obr;
                double y_end = y_beg + h_y_obr;
                Am_obr[I][J] = Integr_Revers(x_beg, x_end, y_beg, y_end, point_view[I][0], point_view[I][1], K0);
            }
            Vec_obr[I] = create_noise(getIntensivity(point_view[I][0], point_view[I][1], A, C, h_x, h_y, n, Vec), .01)
                - fallWave(K0, point_view[I][0], point_view[I][1], x_ist, y_ist);
        }

        //cout << "Cond = " << cond(Am_obr) << endl;
        //----------преобуславливание---------
        //complex<double>** conj_Am, * conj_Vec;
        //CreateMatrix(N, conj_Am);
        //CreateVec(N, conj_Vec);

        //conjMatrix(Am_obr, conj_Am, N);

        //complex<double>** symmetric_Am, * symmetric_Vec;
        //symmetric_Am = mult_matr_matr(conj_Am, Am_obr, N);
        //symmetric_Vec = mult_matr_vec(conj_Am, Vec_obr, N);

        //
        ////complex<double> tmp = mult_matr_matr()
        //complex<double>** precond = predobusl(symmetric_Am, N, .1);
        //complex<double>** preAm = mult_matr_matr(precond, symmetric_Am, N), * preVec = mult_matr_vec(precond, symmetric_Vec, N);
        //cout << "Cond = " << cond(preAm) << endl;
        ////Gauss(preAm, preVec, N_obr);
        ////Vec_obr = BICGstab(preAm, preVec, N);
        //complex<double>* W;
        //CreateVec(N, W);
        //for (size_t i = 0; i < N; i++)
        //{
        //    W[i] = 1.0;
        //}
        //Gradient(N, Am_obr, Vec_obr, W, Vec_obr);
        //----------преобуславливание---------
        Gauss(Am_obr, Vec_obr, N_obr);
        //printInFile(Vec_obr, A, C, h_x_obr, h_y_obr, n_obr, "vosst_alpha.txt");
        complex<double>* vosst_k, * ish_k;
        CreateVec(N_obr, vosst_k);
        CreateVec(N_obr, ish_k);

        for (size_t I = 0; I < N_obr; I++)
        {
            int i_koll = I / n_obr;
            int j_koll = I % n_obr;
            double x_koll = A + i_koll * h_x_obr + h_x_obr / 2.;
            double y_koll = C + j_koll * h_y_obr + h_y_obr / 2.;

            complex<double> Int = fallWave(K0, x_koll, y_koll, x_ist, y_ist);
            for (size_t J = 0; J < N_obr; J++)
            {
                int i_int = J / n_obr;
                int j_int = J % n_obr;
                double x_beg = A + i_int * h_x_obr;
                double y_beg = C + j_int * h_y_obr;
                double x_end = x_beg + h_x_obr;
                double y_end = y_beg + h_y_obr;
                Int += Integr_Revers(x_beg, x_end, y_beg, y_end, x_koll, y_koll, K0) * Vec_obr[J];
            }
            complex<double> tmp = Vec_obr[I] / Int + K0 * K0;
            vosst_k[I] = tmp;
            ish_k[I] = pow(K_f(i_koll, j_koll), 2);
        }
        main_file << "[";
        for (size_t i = 0; i < N_obr; i++)
        {
            main_file << vosst_k[i].real();
            if (i != N_obr - 1)
                main_file << ", ";
        }
        main_file << "]\n";

        //printInFile(vosst_k, A, C, h_x_obr, h_y_obr, n_obr, "vosst_k"+to_string(iteration)+".txt");
        //printInFile(filter(vosst_k, N, iteration+1), A, C, h_x_obr, h_y_obr, n_obr, "vosst_k.txt");
        //printInFile(vosst_k, A, C, h_x_obr, h_y_obr, n_obr, "vosst_k.txt");
        printInFile(ish_k, A, C, h_x_obr, h_y_obr, n_obr, "ish_k.txt");

        ans << "[";
        for (size_t i = 0; i < N_obr; i++)
        {
            ans << ish_k[i].real();
            if (i != N_obr - 1)
                ans << ", ";
        }
        ans << "]\n";
        //return 0;
    }
    ans.close();
    main_file.close();
    MPI_Finalize();
    return 0;
}