#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using Eigen::MatrixXd;
using namespace std;

int main() {
    // 初始化矩阵
    MatrixXd m(2, 2);
    m(0, 0) = 3;
    m(1, 0) = 2.5;
    m(0, 1) = -1;
    m(1, 1) = m(1, 0) + m(0, 1);
    cout << "Matrix m:\n" << m << endl;

    // 外部数据数组
    double mat2[4] = {1.0, 2.0, 3.0, 4.0};

    // 使用Eigen::Map将一维数组转换为矩阵
    MatrixXd test2 = Eigen::Map<MatrixXd>(mat2, 2, 2);
    cout << "Matrix test2:\n" << test2 << endl;

    // 使用SelfAdjointEigenSolver计算特征值和特征向量
    Eigen::SelfAdjointEigenSolver<MatrixXd> eigensolver(test2);
    if (eigensolver.info() != Eigen::Success) {
        cerr << "计算特征值和特征向量失败!" << endl;
        return -1;
    }

    // 输出特征值
    cout << "特征值:\n" << eigensolver.eigenvalues() << endl;

    // 输出特征向量
    cout << "特征向量:\n" << eigensolver.eigenvectors() << endl;

    return 0;
}