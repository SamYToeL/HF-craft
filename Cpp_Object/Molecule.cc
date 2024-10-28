#include "Molecule.h"
#include<iostream>
#include<fstream>
#include <cstdio>
using namespace std;



double squa(double x){
    return x*x;
}

Molecule::Molecule(const char *filename, int q){
    int charge = q;
    ifstream tar(filename);
    assert(tar.good());

    //Set natom
    tar >> natom;
    charge = q;
    
    //allocate space
    ele = new int[natom];
    geom = new double *[natom];      //地址中的地址 double **geom    geom = new double *[ele]

    for(int i=0;i<natom;i++){
        geom[i] = new double[3];
        tar >> ele[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
    }
    tar.close();

}

void Molecule::print_geom(){
    
    cout << natom << endl;
    //printf("%d\n",natom);
    for(int i=0;i<natom; i++){
        printf("%d %8.5f %8.5f %8.5f\n", ele[i],geom[i][0],geom[i][1],geom[i][2]);
    }
}

double Molecule::bond(int atom1, int atom2){
    double distance;
    distance = sqrt(squa(geom[atom1][0]-geom[atom2][0])+squa(geom[atom1][1]-geom[atom2][1])+squa(geom[atom1][3]-geom[atom2][3]));
    return distance;
}

Molecule::~Molecule(){
    delete[] ele;
    for(int i=0;i<natom;i++){
        delete[] geom[i];
    }
    delete geom;
}