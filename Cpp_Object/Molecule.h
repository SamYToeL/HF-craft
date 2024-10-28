#include<string>
using namespace std;

class Molecule{
    public:
        int natom;
        int charge;
        int *ele;
        double **geom;
        void print_geom();
        double bond(int atom1, int atom2);
        double angle(int atom1, int atom2, int atom3);
        double torsion(int atom1, int atom2, int atom3, int atom4);
    
        Molecule(const char *filename, int q);
        ~Molecule();
};


