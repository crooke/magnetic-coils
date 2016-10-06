/* MagneticField.cpp
 * Author: Eric N. Crook
 * Date: April 24, 2016
 * Description: Implementation of the MagneticField class.
 * Used to determine the magnetic field produced by a coil(s) of wire.
*/

#include "header.h"
#include <boost/math/special_functions/ellint_1.hpp> //boost package for elliptic integrals
#include <boost/math/special_functions/ellint_2.hpp>

MagneticField::MagneticField(double _rMax, int _nr, double _zMin, double _zMax, int _nz) {
    //assign function parameters to private variables
    rMax = _rMax;
    zMin = _zMin;
    zMax = _zMax;
    nr = _nr;
    nz = _nz;

    //create two arrays to hold r and z coordinates/points
    r = new double[nr];
    z = new double[nz];

    //fill the array with the points based on the given number of points
    // we want in the given range
    double dr = rMax / (nr-1); //distance b/w r points
    double dz = (zMax - zMin) / (nz-1);

    for (int i = 0; i < nr; i++)
        r[i] = i * dr;

    for (int j = 0; j < nz; j++)
        z[j] = j * dz + zMin;

    //create the Br and Bz arrays
    Br = new double* [nr];

    for (int i = 0; i < nr; i++) {
        Br[i] = new double[nz];
    }

    Bz = new double* [nr];

    for (int i = 0; i < nr; i++)
        Bz[i] = new double[nz];

    //initialize Br and Bz elements to zero
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nz; j++) {
            Br[i][j] = 0.0;
            Bz[i][j] = 0.0;
        }
    }

}

MagneticField::~MagneticField() {
    delete [] r; //destroys (freeing up heap/freestore) array r was pointing to
    r = nullptr; //r points to null pointer

    delete [] z;
    z = nullptr;

    for (int i = 0; i < nr; i++) {
            delete [] Br[i]; //delete the array that the ith element of Br points to
            delete [] Bz[i];
    }

    delete [] Br;
    delete [] Bz;

    Br = nullptr;
    Bz = nullptr;
}

void MagneticField::coils(double _innerR, double _outerR, double _width, double _dist, double _current, double _gauge, double _separation) {
    //Define the Coil
    //assign parameters to private members:
    innerR = _innerR;
    outerR = _outerR;
    width = _width;
    dist = _dist;
    I = _current;
    gauge = _gauge;
    separation = _separation;
}

void MagneticField::calcOnePt(double rr, double zz) {
    double d = gauge + separation; //distance b/w the wire loops
    int nlr = int ((outerR - innerR)/d); //number of loops in the coil (radial direction)
    int nlz = int (width/d); //number of loops in the coil's width (z direction)
    double a; //radius loops

    int ri = rToIndex(rr);
    int zi = zToIndex(zz);

    for (int p = 0; p < nlr; p++) { //sum over all the loops in the radial direction
        a = innerR + p * d; //radius of successive loops

        for (int q = 0; q < nlz+1; q++) { //sum over all the loops in the z-direction
            if (a != r[ri]) { //prevents elliptic function from throwing error when k = 1
                //top coil
                Br[ri][zi] += calcBr(rr, zz - dist/2 - q * d, a);
                Bz[ri][zi] += calcBz(rr, zz - dist/2 - q * d, a);
                //bottom coil
                Br[ri][zi] += calcBr(rr, zz + dist/2 + q * d, a);
                Bz[ri][zi] += calcBz(rr, zz + dist/2 + q * d, a);
            }
        }
    }

    cout << Br[ri][zi] << "\t" << Bz[ri][zi] << endl;
}

void MagneticField::calculate() {
    double d = gauge + separation; //distance b/w the wire loops
    int nlr = int ((outerR - innerR)/d); //number of loops in the coil (radial direction)
    int nlz = int (width/d); //number of loops in the coil's width (z direction)

    double h = dist / 2; //distance from origin (halfway between top and bottom coil) 
                         // to the first current loop
    double a; //radius of loops
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nz; j++) {
            for (int p = 0; p < nlr; p++) { //sum over all the loops in the radial direction
                a = innerR + p * d; //radius of successive loops

                for (int q = 0; q < nlz; q++) { //sum over all the loops in the z-direction
                    if (a != r[i]) { //prevents elliptic func. from throwing error when k = 1
                        //top coil
                        Br[i][j] += calcBr(r[i], z[j] - h - q * d, a);
                        Bz[i][j] += calcBz(r[i], z[j] - h - q * d, a);
                        //bottom coil
                        Br[i][j] += calcBr(r[i], z[j] + h + q * d, a);
                        Bz[i][j] += calcBz(r[i], z[j] + h + q * d, a);
                    }
                }
            }
        }
    }
}

void MagneticField::outputField(string file_name) {
    ofstream out;
    out.open(file_name);

    out << "#r\tz\tBr\tBz\n";

    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nz; j++) {
            //print B where the coil is NOT located
            if (!((i > rToIndex(innerR) && i < rToIndex(outerR) && j > zToIndex(dist/2) && j < zToIndex(dist/2 + width))
                  || (i > rToIndex(innerR) && i < rToIndex(outerR) && j < zToIndex(-dist/2) && j > zToIndex(-dist/2 - width)) )) {
                out << r[i] << "\t" << z[j] << "\t" << Br[i][j] << "\t" << Bz[i][j] << 0 << endl;
            }
            else
                out << r[i] << "\t" << z[j] << "\t" << 0 << "\t" << 0 << "\t" << endl;
        }
    }
    out.close();
}

void MagneticField::outputGeometry(string file_name) {
    ofstream out;
    out.open(file_name);

    out << "#r\tz\tcolor (1 - coils, 0 - elsewhere)\n";

    //geometry of coil; outputs so gnuplot can color it in
    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nz; j++) {
            if (!(i > rToIndex(innerR) && i < rToIndex(outerR) && j > zToIndex(dist/2) && j < zToIndex(dist/2 + width))) {
                out << r[i] << "\t" << z[j] << "\t" << 0 << endl;
            }
            else
                out << r[i] << "\t" << z[j] << "\t" << 1 << endl;
        }
    }
    out.close();
}

double MagneticField::calcBr(double r, double z, double a) {
    using namespace boost::math;
    double k_2 = (4 * r * a) / (z * z + (a + r) * (a + r)); //k^2
    return (MU * I * z)/(2 * PI * r * sqrt(z*z + (a + r)*(a + r))) * ((a*a + z*z + r*r)/(z*z + (r - a)*(r - a)) * ellint_2(k_2) - ellint_1(k_2));
}

double MagneticField::calcBz(double r, double z, double a) {
    using namespace boost::math;
    double k_2 = (4 * r * a) / (z * z + (a + r) * (a + r)); //k^2
    return (MU * I)/(2 * PI * sqrt(z*z + (a + r)*(a + r))) * ((a*a - z*z - r*r)/(z*z + (r - a)*(r - a)) * ellint_2(k_2) + ellint_1(k_2));
}

int MagneticField::rToIndex(double rr) {
    return coordToIndex(rr, 0, rMax, nr);
}

int MagneticField::zToIndex(double zz) {
    return coordToIndex(zz, zMin, zMax, nz);
}

int MagneticField::coordToIndex(double cc, double cMin, double cMax, int nc) {
    int i = int ((cc - cMin)*(nc - 1)/(cMax - cMin));
    return i; }
