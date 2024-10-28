#include<iostream>
#include<string>
#include "Molecule.h"
using namespace std;

int main(){
    Molecule Tar_mol("allene.dat",0);
    int num = Tar_mol.natom;
    for(int i = 0;i < num; i++){
        for(int j = i+1; j < num; j++){
            float d = Tar_mol.bond(i,j);
            cout << "the distance between atom: " << i+1 << " and atom " << j+1 << " is " << d << endl;
        }
    }
    return 0;
}
