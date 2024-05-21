#pragma once
#include <cstring>

#ifdef _NEC
#include <aslfftw3.h>
#else
#include <fftw3.h>
#endif

class FFTW_Executor {
private:
    int Nx = 0;
    int Ny = 0;
    int Nz = 0;    
    fftw_complex* buffer = nullptr;
    fftw_plan plan_forward{ 0 };
    fftw_plan plan_backward{ 0 };
public:
    ~FFTW_Executor() {
        if (buffer) {
            fftw_destroy_plan(plan_forward);
            fftw_destroy_plan(plan_backward);
            fftw_free(buffer);
        }
    }

    void Initialize(int grid_x, int grid_y, int grid_z, int fftw_flags) {
        Nx = grid_x;
        Ny = grid_y;
        Nz = grid_z;
        buffer = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx * Ny * Nz);
        plan_forward = fftw_plan_dft_3d(Nz, Ny, Nx, buffer, buffer, FFTW_FORWARD, fftw_flags);
        plan_backward = fftw_plan_dft_3d(Nz, Ny, Nx, buffer, buffer, FFTW_BACKWARD, fftw_flags);
    }

    fftw_complex* GetBuffer() { return buffer; };

    fftw_complex* ForwardExecute(fftw_complex* in, fftw_complex* out) {
        if (in && (in != buffer)) {
            std::memcpy(buffer, in, sizeof(fftw_complex*) * Nx * Ny * Nz);
        }
        fftw_execute(plan_forward);
        if (out && (out != buffer)) {
            std::memcpy(out, buffer, sizeof(fftw_complex*) * Nx * Ny * Nz);
        }
        return buffer;
    }

    fftw_complex* BackwardExecute(fftw_complex* in, fftw_complex* out) {
        if (in && (in != buffer)) {
            std::memcpy(buffer, in, sizeof(fftw_complex*) * Nx * Ny * Nz);
        }
        fftw_execute(plan_backward);
        if (out && (out != buffer)) {
            std::memcpy(out, buffer, sizeof(fftw_complex*) * Nx * Ny * Nz);
        }
        return buffer;
    }

};
