#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <Eigen/Dense>

using namespace std;

float *read_NNrep(const char *file) {
    ifstream rep(file, ios::in);
    if (!rep.is_open()) {
        cerr << "Unable to open file " << file << endl;
        return nullptr;
    }
    float *reps = new float[100];
    int i = 0;
    while (rep >> reps[i] && i < 100) {
        i++;
    }
    return reps;
}

double **read_Smat(const char *file) {
    ifstream Sfile(file, ios::in);
    if (!Sfile.is_open()) {
        cerr << "Unable to open file " << file << endl;
        return nullptr;
    }
    double **S_inte = new double *[7];
    for (int k = 0; k < 7; k++) {
        S_inte[k] = new double[7]();
    }
    while(!Sfile.eof()) {
        int atom1, atom2;
        double value;
        Sfile >> atom1 >> atom2 >> value;
        S_inte[atom1-1][atom2-1] = value;
        S_inte[atom2-1][atom1-1] = value;
    }
    cout << S_inte << endl;
    return S_inte;
}

float **read_Tmat(const char *file) {
    ifstream Tfile(file, ios::in);
    if (!Tfile.is_open()) {
        cerr << "Unable to open file " << file << endl;
        return nullptr;
    }
    float **T_inte = new float *[7];
    for (int k = 0; k < 7; k++) {
        T_inte[k] = new float[7]();
    }
    for (int i = 0; i < 28; i++) {
        int atom1, atom2;
        float value;
        Tfile >> atom1 >> atom2 >> value;
        if (Tfile.fail()) {
            cerr << "Error reading file " << file << " at line " << i + 1 << endl;
            break;
        }
        T_inte[atom1][atom2] = value;
        T_inte[atom2][atom1] = value;
    }
    return T_inte;
}

float **read_vmat(const char *file) {
    ifstream vfile(file, ios::in);
    if (!vfile.is_open()) {
        cerr << "Unable to open file " << file << endl;
        return nullptr;
    }
    float **v_inte = new float *[7];
    for (int k = 0; k < 7; k++) {
        v_inte[k] = new float[7]();
    }
    for (int i = 0; i < 28; i++) {
        int atom1, atom2;
        float value;
        vfile >> atom1 >> atom2 >> value;
        if (vfile.fail()) {
            cerr << "Error reading file " << file << " at line " << i + 1 << endl;
            break;
        }
        v_inte[atom1][atom2] = value;
        v_inte[atom2][atom1] = value;
    }
    return v_inte;
}

int *find_mat(int num) {
    int *ofw = new int[num];
    ofw[0] = 0;
    for (int i = 1; i < num; i++) {
        ofw[i] = ofw[i - 1] + i;
    }
    return ofw;
}

float *read_twoint(const char *file) {
    ifstream ingfile(file, ios::in);
    if (!ingfile.is_open()) {
        cerr << "Unable to open file " << file << endl;
        return nullptr;
    }
    int i, j, p, q, ij, pq, ijpq;
    float value;
    int max_size = 665;
    float *trray = new float[max_size]();
    int *ioff = find_mat(1000);

    while (ingfile >> i >> j >> p >> q >> value) {
        ij = (i > j) ? ioff[i] + j : ioff[j] + i;
        pq = (p > q) ? ioff[p] + q : ioff[q] + p;
        ijpq = (ij > pq) ? ioff[ij] + pq : ioff[pq] + ij;
        trray[ijpq] = value;
    }
    delete[] ioff;
    return trray;
}

Eigen::MatrixXd get_inver_sqrt(Eigen::MatrixXd mat, int size) {
    Eigen::MatrixXd inv = mat.inverse();
    for (int i = 0; i < size; i++) {
        inv(i, i) = sqrt(inv(i, i));
    }
    return inv;
}

int main() {
    double **S_mat = read_Smat("s.dat");
    if (S_mat == nullptr) {
        return -1;
    }
    
    printf("so far so good");
    Eigen::MatrixXd S_am(7, 7);
    for (int i = 0; i < 7; i++) {
        for (int j = 0; j < 7; j++) {
            S_am(i, j) = S_mat[i][j];
        }
    }

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S_am);
    if (eigensolver.info() != Eigen::Success) {
        cerr << "Eigenvalue computation failed!" << endl;
        return -1;
    }

    Eigen::MatrixXd egmt = eigensolver.eigenvectors();
    Eigen::MatrixXd egvs = eigensolver.eigenvalues().asDiagonal();
    Eigen::MatrixXd egv = get_inver_sqrt(egvs, 7);
    Eigen::MatrixXd S_minus5 = egmt * egv * egmt.transpose();

    cout << "eigen value matrix is \n" << egvs << endl;
    
    //cout << S_minus5 << endl;
    
    for (int i = 0; i < 7; i++) {
        delete[] S_mat[i];
    }
    delete[] S_mat;

    return 0;
}
