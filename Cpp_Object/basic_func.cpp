#include<iostream>
#include "Molecule.cc"
#include<string>
using namespace std;

int main(){
    const char *pt = "/Users/samliao/Desktop/Code/Cpp_Object/Object1/benzene.dat";
    Molecule testm(pt,0);
    testm.print_geom();
    return 0;   
}