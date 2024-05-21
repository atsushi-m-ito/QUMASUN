#pragma once
#include "GridRange.h"
#include "GridFor.h"
#include "vec3.h"
#include "soacomplex.h"

inline
void AddVelocityForWave(const GridRange& grid, SoAComplex&& psi, vec3d v, vec3d r0, const GridRange& global_grid, double dx, double dy, double dz){
    
    const double box_x = global_grid.SizeX() * dx;
    const double box_y = global_grid.SizeY() * dy;
    const double box_z = global_grid.SizeZ() * dz;
    const double half_x = box_x / 2.0;
    const double half_y = box_y / 2.0;
    const double half_z = box_z / 2.0;


    ForXYZ(grid, [&](int64_t i, int64_t ix, int64_t iy, int64_t iz) {
        double x = (double)ix * dx - r0.x;
        double y = (double)iy * dy - r0.y;
        double z = (double)iz * dz - r0.z;
        if (x > half_x) x -= box_x;
        if (x < -half_x) x += box_x;
        if (y > half_y) y -= box_y;
        if (y < -half_y) y += box_y;
        if (z > half_x) z -= box_z;
        if (z < -half_x) z += box_z;

        double vx = v.x * x + v.y * y + v.z * z;
        double psi0_r = psi.re[i];
        double psi0_i = psi.im[i];
        const double cos_vx = cos(vx);
        const double sin_vx = sin(vx);
        psi.re[i] = cos_vx * psi0_r - sin_vx * psi0_i;
        psi.im[i] = cos_vx * psi0_i + sin_vx * psi0_r;

        });
}
