#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

#define INDEX(i,j) (i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i)

double read_NNrep(const char *file){ 
    ifstream rep(file, ios::in);
    double reps;
    rep >> reps;
    return reps;
    //记得在main最后删除该数组
}

double **read_Smat(const char *file){
    ifstream Sfile(file, ios::in);
    double **S_inte = new double *[7];
    for(int k=0;k<7;k++){
        S_inte[k] = new double [7];
    }
    for(int i = 0;i< 28;i++){
        int atom1, atom2;
        double value;
        Sfile >> atom1 >> atom2 >> value;
        if(atom1 != atom2){
            S_inte[atom2-1][atom1-1] = value;
            S_inte[atom1-1][atom2-1] = value;
        }
        else{
            S_inte[atom1-1][atom2-1] = value;
        }//not a necessary step to do this if...else judgement.
    }

    return S_inte;
}  
    //记得loop删除该指针，and 取值的时候顺序要记得调换
    

double **read_Tmat(const char *file){
    ifstream Tfile(file, ios::in);
    double **T_inte = new double*[7];
    for(int k=0;k<7;k++){
        T_inte[k] = new double [7];
    }
    for(int i = 0;i< 28;i++){
        
        int atom1, atom2;
        double value;
        Tfile >> atom1 >> atom2 >> value;
        T_inte[atom1-1][atom2-1] = value;
        T_inte[atom2-1][atom1-1] = value;
    } 
    //记得loop删除该指针，and 取值的时候顺序要记得调换 index 大的在前面
    return T_inte;
}

double **read_Vmat(const char *file){
    ifstream vfile(file, ios::in);
    double **v_inte = new double *[7];
    for(int k=0;k<7;k++){
        v_inte[k] = new double [7];
    }    
    for(int i = 0;i< 28;i++){
        int atom1, atom2;
        double value;
        vfile >> atom1 >> atom2 >> value;
        v_inte[atom1-1][atom2-1] = value;
        v_inte[atom2-1][atom1-1] = value;
    }
    //记得loop删除该指针，and 取值的时候顺序要记得调换 index 大的在前面
    return v_inte;
}

int *find_mat(int num){
    int *ofw = new int[num];
    ofw[0] = 0;
    for(int i=1;i<num;i++)
        ofw[i] = ofw[i-1] + i;
    return ofw;
    //Define a locating matrix by recursion
}


int locn(int i,int j,int p,int q,int* ioff){
    int ij,pq,ijpq;
    ij = (i>j) ? (ioff[i] + j): (ioff[j] + i);
    pq = (p>q) ? (ioff[p] + q): (ioff[q] + p);
    ijpq = (ij > pq) ? (ioff[ij] + pq): (ioff[pq] + ij);
    return ijpq;
}

double *read_twoint(const char *file){
    ifstream ingfile(file, ios::in);
    int i,j,k,l,ij,kl,ijkl;
    double value;
    int max_size = 1000;
    double *trray = new double[max_size];
    int *ioff = find_mat(1000);

    while(!ingfile.eof()){
        ingfile >> i >> j >> k >> l >> value;
        ij = INDEX(i,j);
        kl = INDEX(k,l);
        ijkl = INDEX(ij,kl);
        //cout<< ijkl << endl;
        trray[ijkl] = value;
        cout << i << j << k << l << "   " << ij << ' ' << kl << "  "<< value << endl;
    }
    //cout << count <<endl;
    delete[] ioff;
    return trray;
}

Eigen::MatrixXd get_inver_sqrt(Eigen::MatrixXd mat, int size){
    Eigen::MatrixXd inv = mat.inverse();
    for(int i=0; i < size; i++){
        inv(i,i) = sqrt(inv(i,i));
    }
    return inv;
}

int get_orb(const char *file){
    ifstream tar(file, ios::in);
    int max_orb=0;
    int orb1, orb2;
    double value;
    while(tar>> orb1 >> orb2 >> value){
        if(orb1 > max_orb) max_orb = orb1;
        if(orb2 > max_orb) max_orb = orb2;
    }
    return max_orb;
}

MatrixXd sym_ortho(MatrixXd S_am,int norb){

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(S_am);
    Eigen::MatrixXd egmt = eigensolver.eigenvectors();
    Eigen::MatrixXd egvs = eigensolver.eigenvalues().asDiagonal();
    //Eigen ^-1/2
    Eigen::MatrixXd egv = get_inver_sqrt(egvs, norb);
    //S^-1/2
    Eigen::MatrixXd S_minus5 = egmt*egv*egmt.transpose(); // eigenvectors * eigenvalues ^-0.5 * eigenvectors^Transpoes 
    //cout << "The S^-1/2 is \n" << S_minus5 << endl;
    return S_minus5;

}

MatrixXd Dens_Mat(MatrixXd S_minus5, MatrixXd F_ini,int norb,int elec){
    Eigen::MatrixXd F_0 = S_minus5.transpose()*(F_ini)*S_minus5;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es_F(F_0);
    Eigen::MatrixXd C_p = es_F.eigenvectors();
    Eigen::MatrixXd C_0 = S_minus5*C_p;
    Eigen::MatrixXd D_0(norb,norb);
    D_0 = Eigen::MatrixXd::Zero(norb,norb);  //Not sure if it works
    for(int i=0;i<norb;i++){
        for(int j=0;j<norb;j++){
            for(int m=0;m< elec/2 ;m++){
                D_0(i,j) += C_0(i,m)*C_0(j,m);
            }
            
        }
    }
    //cout << D_0 <<endl;
    return D_0;
}

MatrixXd new_Fock(MatrixXd &H_core, double *trray, MatrixXd &D, int * ioff,int norb){
    //cout << H_core << endl;
    int ij,kl,ijkl,ik,jl,ikjl;
    MatrixXd Fock_m(norb,norb);
    Fock_m = Eigen::MatrixXd::Zero(norb,norb);
    for(int i =0; i< norb;i++)
        for(int j=0; j< norb; j++){
            Fock_m(i,j) += H_core(i,j);
            //double tmp_val = 0.0;
            for(int k =0; k< norb;k++)
                for(int l =0; l< norb;l++){
                    //int loc = locn(m,n,i,j,ioff);
                    //int loc2 = locn(m,i,n,j,ioff);
                    ij = INDEX(i+1,j+1);
                    kl = INDEX(k+1,l+1);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i+1,k+1);
                    jl = INDEX(j+1,l+1);
                    ikjl = INDEX(ik,jl);
                    Fock_m(i,j) += D(k,l)*(2.0 * trray[ijkl] - trray[ikjl]);
                }
        }
    //cout << Fock_m << endl;
    return Fock_m;
}



double SCF_E(MatrixXd D_0, MatrixXd &H_core, MatrixXd F_ini,int norb){
    double E_0 = 0;
    for(int m =0; m< norb;m++){
        for(int n=0; n< norb; n++){
            E_0 += D_0(m,n)*(H_core(m,n)+F_ini(m,n));
        }
        //cout<< "Fock(m,n) after\n" << Fock_m << endl;
    }
    //cout<< "Fockm = \n" << Fock_m << "\n"<< endl;
    //cout<< "Hcore = \n" << H_core << "\n"<< endl;
    return E_0;
}



int main(){
    int elec = 10;
    double rps = read_NNrep("enuc.dat");
    double **S_mat = read_Smat("s.dat");
    double **T_mat = read_Tmat("t.dat");
    double **V_mat = read_Vmat("v.dat");
    double *trray = read_twoint("eri-copy.dat");
    int norb = get_orb("s.dat");
    int *iff = find_mat(1000);
    cout << *trray << endl;
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
    //E_0 = SCF_E(D_0,H_core,F_ini,norb);  //get E_0
    tmp_Fock = new_Fock(H_core, trray, D_0,iff,norb);
    //cout << tmp_Fock << endl;
    //E_total = E_0 + rps;
    int i,j,k,l,ij,kl,ik,jl,ijkl,ikjl;
    i=1;
    j=0;
    k=0;
    l=0;
    int i1 = i+1;
    int j1 = j+1;
    int k1 = k+1;
    int l1 = l+1;
                    ij = INDEX(i1,j1);
                    kl = INDEX(k1,l1);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i1,k1);
                    jl = INDEX(j1,l1);
                    ikjl = INDEX(ik,jl);
                    //Fock_m(i,j) += D(k,l)*(2.0 * trray[ijkl] - trray[ikjl]);
    cout << ij << kl << "   "<< ijkl << "  : " << trray[ijkl] << endl;
    cout << H_core(0,0)<< endl;

    return 0;
}