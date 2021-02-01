#ifndef thermocouplesPlane_H
#define thermocouplesPlane_H
#define _USE_MATH_DEFINES

#include "ITHACAutilities.H"

/// Class defining a plane for thermocouples measurements interpolation
class thermocouplesPlane
{
    public:
        double minX;
        double minZ;
        double Y;
        double maxX;
        double maxZ;
        List<scalar> Tmeas;
        List<scalar> Tdirect;
        List<scalar> Tsens;
        List<scalar> Tdiff;
        List<scalar> cellVol;
        List<int> cellID;
        List<scalar> thermocoupleX;
        List<scalar> thermocoupleZ;
        Foam::vector thermocoupleCellDim;
};
#endif
