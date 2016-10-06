/* MagneticField.h
 * Author: Eric N. Crook
 * Date: April 24, 2016
 * Description: Interface for the MagneticField class.
 * Used to determine the magnetic field produced by a coil(s) of wires.
*/

#ifndef RESEARCHCODE_HEADER_H
#define RESEARCHCODE_HEADER_H

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

class MagneticField {
public:
    MagneticField(double _rMax, int _nr, double _zMin, double _zMax, int _nz);
    //creates the 2D space (since coil of wire is symmetric, no azimuthal component)
    // we want to know mag field over

    ~MagneticField();
    //destructor

    void coils(double _innerR, double _outerR, double _width, double _dist, double _current, double _gauge, double _separation);
    //defines the coils being modeled

    void calcOnePt(double rr, double zz);

    void calculate();
    //calculate magnetic field

    void outputField(string file_name);
    //outputs magnetic field values into a file that can be plotted

    void outputGeometry(string file_name);
    //outputs values for the geometry of the coils into a file that can be plotted

private:
    //functions:
    double calcBr(double r, double z, double a);
    //input a point (r,z) and radius a of current loop;
    // returns r-comp of mag field at that point
    //units: input--meters, output--Tesla

    double calcBz(double r, double z, double a);
    //input a point (r,z) and radius a of current loop;
    //returns z-comp of mag field at that point

    int rToIndex(double rr); //input a point r (rr) and returns the index it should be assigned to
    int zToIndex(double zz);
    int coordToIndex(double cc, double cMin, double cMax, int nc); //returns the index a point
    //cc should be assigned to based on cMin, etc.

    //variables:
    double rMax; //maximum radial distance away from the origin we want to calculate B
    //for(cylindrical coordinates)
    int nr; //number of r points
    double zMin; //minimum z
    double zMax; //maximum z
    int nz; //number of z points

    double *r; //points to an array that holds the r points
    double *z; //points to an array that holds the z points

    double **Br; //points to a two-dim array that holds the r-comp of mag field at points (r,z)
    double **Bz; // ^                                       z-comp


    double innerR; //inner radius of coils
    double outerR; //outer radius of coils
    double width; //width or thickness of the magnets
    double dist; //distance b/w coils
    double I; //current
    double gauge; //diameter of wires in meters (not actually AWG)
    double separation; //distance b/w wires in the coil

    const double PI = 3.1459;
    const double MU = 4 * PI * pow(10,-7);
};

#endif //RESEARCHCODE_HEADER_H
