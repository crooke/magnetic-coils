/* main.cpp
 * Author: Eric N. Crook
 * Date: April 24, 2016
 * Description: This application uses the MagneticField class in order to calculate 
 * and ouptut the magnetic field produced by the coils currently used in the 
 * negative ion matter interferometry experiment.
*/

#include <iostream>
#include "MagneticField.h"

using namespace std;

int main() {
    //Actual lab coil:  (don't forget, units are in meters!)
    MagneticField LabMagnets = MagneticField(.4, 21, -0.2, 0.2, 21);
    LabMagnets.coils(0.1113, 0.2483, 0.0635, 0.1349, -1, 0.0020525, 0);
    LabMagnets.calculate();
    LabMagnets.outputField("ResearchCodeMagneticField-ZoomedOut.data");
    LabMagnets.outputGeometry("geometry.data");

    //"zoomed in" view of the drift region
    MagneticField DriftRegion = MagneticField(.4, 41, -0.2, 0.2, 41);
    DriftRegion.coils(0.1113, 0.2483, 0.0635, 0.1349, -1, 0.0020525, 0);
    DriftRegion.calculate();
    DriftRegion.outputField("ResearchCodeMagneticField-Zoomed-In.data");

    //used compare the measurements taken with the gaussmeter to the programs calculations
    DriftRegion.calcOnePt(0.0254, 0);
    DriftRegion.calcOnePt(0.0508, 0);
    DriftRegion.calcOnePt(0.0762, 0);
    DriftRegion.calcOnePt(0.1016, 0);
    DriftRegion.calcOnePt(0.1270, 0);
    DriftRegion.calcOnePt(0.1524, 0);
    DriftRegion.calcOnePt(0.1778, 0);
    DriftRegion.calcOnePt(0.2032, 0);
    DriftRegion.calcOnePt(0.2286, 0);
    DriftRegion.calcOnePt(0.2540, 0);

    return 0;
}
