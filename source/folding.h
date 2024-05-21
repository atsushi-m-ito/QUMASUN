#pragma once
#include <cmath>

inline 
void Folding(double& rx, double& ry, double& rz, const double* box_axis) {
    double ix = rx / box_axis[0];
    ix -= std::floor(ix);
    rx = ix * box_axis[0];

    double iy = ry / box_axis[4];
    iy -= std::floor(iy);
    ry = iy * box_axis[4];

    double iz = rz / box_axis[8];
    iz -= std::floor(iz);
    rz = iz * box_axis[8];
}
