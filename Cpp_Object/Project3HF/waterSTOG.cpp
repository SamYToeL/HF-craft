//#include<iostream>
//#include<fstream>
//#include<cmath>
//#include<string>
//#include<Eigen/Dense>
#include"Functions.cc"
using namespace std;
using namespace Eigen; 
//This project followed tutorials updated by Prof. D.Crawford. at  
//https://github.com/CrawfordGroup/ProgrammingProjects/tree/master

int main(){
    int elec = 10;
    double rps = read_NNrep("enuc.dat");
    double **S_mat = read_Smat("s.dat");
    double **T_mat = read_Tmat("t.dat");
    double **V_mat = read_Vmat("v.dat");
    double *trray = read_twoint("eri.dat");
    int norb = get_orb("s.dat");
    int *iff = find_mat(1000);
    //补丁 at first I store the integral in arrays, which is not convenient to operate,
    //but so fat I dont want to change those functions, thus the next coming are three patches
    Eigen::MatrixXd S_am(norb,norb);
    for(int i = 0;i< norb;i++){
        for(int j = 0;j< norb;j++){
            S_am(i,j) = S_mat[i][j];
        }
    }
    Eigen::MatrixXd T_m(norb,norb);
    for(int i = 0;i< norb;i++){
        for(int j = 0;j< norb;j++){
            T_m(i,j) = T_mat[i][j];
        }
    }
    Eigen::MatrixXd V_m(norb,norb);
    for(int i = 0;i< norb;i++){
        for(int j = 0;j< norb;j++){
            V_m(i,j) = V_mat[i][j];
        }
    }

    MatrixXd S_minus5 = sym_ortho(S_am,norb);

    //Then give the initial guess of density matrix
    Eigen::MatrixXd H_core = V_m + T_m;
    Eigen::MatrixXd F_ini = H_core;   //initial guess
    MatrixXd D_0 = Dens_Mat(S_minus5,F_ini,norb,elec); //MatrxXd S_minus5, MatrxXd F_ini,int norb,int elec
    //Compute energy & new_fock Sum_m (P_uv*(H^core_uv + F_uv)) // Closed-shell hartree fock   2 and 1/2 cancel out by the way we define density matrix
    MatrixXd tmp_Fock;
    double E_0;
    double E_old, E_tmp;
    double del_E;
    double rmsd;
    MatrixXd D_tmp;
    MatrixXd D_old;
    int round = 1;
    double E_total;
    E_0 = SCF_E(D_0,H_core,F_ini,norb);  //get E_0
    tmp_Fock = new_Fock(H_core, trray, D_0,iff,norb);
    //cout << tmp_Fock << endl;
    E_total = E_0 + rps;
    printf("Round        E(elec)              E(tot)               Delta(E)             RMS(D)\n");    
    cout << round << "   " << "E_0 ="  << E_0  << "   E_0(total) = " << E_total << endl;
    D_old = D_0;
    E_old = E_0;
    //check convergence
    rmsd = 2.0;
    while(del_E > (1e-6) || rmsd > (1e-4)){
        D_tmp = Dens_Mat(S_minus5, tmp_Fock,norb,elec);
        E_tmp = SCF_E(D_tmp,H_core,tmp_Fock,norb);
        E_total = E_tmp + rps;
        rmsd = 0.0;
        for(int i =0; i< norb;i++){
            for(int j =0; j< norb;j++){
                rmsd += (D_tmp(i,j)-D_old(i,j))*(D_tmp(i,j)-D_old(i,j));
            }
        }
        del_E = E_tmp - E_old;
        rmsd = sqrt(rmsd);
        E_old = E_tmp;
        D_old = D_tmp;
        tmp_Fock = new_Fock(H_core, trray, D_tmp,iff,norb);
        round++;
        if (round > 34) {
            cout<< "Conver failure" <<endl;
            break;
        }
           printf("%02d    %.12f    %.12f    %.12f    %.12f\n", 
                    round, E_tmp, E_total, del_E, rmsd);
    }
    cout<< "++++++++++++++++++= iteration end +++++++++++++++++++++++++++++++++" <<endl;
    
    
    for(int i=0; i<7;i++){
        delete[] S_mat[i];  
        delete[] T_mat[i];
        delete[] V_mat[i];
    }
    delete[] S_mat;
    delete[] T_mat;
    delete[] V_mat;
    
    return 0;
}
