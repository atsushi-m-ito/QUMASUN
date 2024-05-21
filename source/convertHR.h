#pragma once
#include <cstdint>

template<class T>
void UpConvert_Kspace(const T* src, int size_x, int size_y, int size_z, T* hr_dest, int hr_ratio_x, int hr_ratio_y, int hr_ratio_z) {

    memset(hr_dest, 0, sizeof(T) * (int64_t)size_x * (int64_t)size_y * (int64_t)size_z * (int64_t)hr_ratio_x * (int64_t)hr_ratio_y * (int64_t)hr_ratio_z);

    for (int iz = 0; iz < size_z/2; ++iz) {
        const int64_t oiz = iz * hr_ratio_x * size_x * hr_ratio_y * size_y;
        for (int iy = 0; iy < size_y / 2; ++iy) {
            const int64_t oiy = iy * hr_ratio_x * size_x;
            for (int ix = 0; ix < size_x / 2; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
            for (int ix = size_x / 2; ix < size_x; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = size_x * (hr_ratio_x - 1) + ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
        }
        for (int iy = size_y / 2; iy < size_y; ++iy) {
            const int64_t oiy = (size_y * (hr_ratio_y - 1) + iy) * hr_ratio_x * size_x;
            for (int ix = 0; ix < size_x / 2; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
            for (int ix = size_x / 2; ix < size_x; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = size_x * (hr_ratio_x - 1) + ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
        }
    }
    for (int iz = size_z / 2; iz < size_z; ++iz) {
        const int64_t oiz = (size_z* (hr_ratio_z-1) + iz) * hr_ratio_x * size_x * hr_ratio_y * size_y;
        for (int iy = 0; iy < size_y / 2; ++iy) {
            const int64_t oiy = iy * hr_ratio_x * size_x;
            for (int ix = 0; ix < size_x / 2; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
            for (int ix = size_x / 2; ix < size_x; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = size_x * (hr_ratio_x - 1) + ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
        }
        for (int iy = size_y / 2; iy < size_y; ++iy) {
            const int64_t oiy = (size_y * (hr_ratio_y - 1) + iy) * hr_ratio_x * size_x;
            for (int ix = 0; ix < size_x / 2; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
            for (int ix = size_x / 2; ix < size_x; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                const int64_t oi = size_x * (hr_ratio_x - 1) + ix + oiy + oiz;
                hr_dest[oi] = src[i];
            }
        }
    }
}


template<class T>
void DownConvert_Realspace(T* dest, int size_x, int size_y, int size_z, const T* hr_src, int hr_ratio_x, int hr_ratio_y, int hr_ratio_z) {
    for (int iz = 0; iz < size_z; ++iz) {
        for (int iy = 0; iy < size_y; ++iy) {
            for (int ix = 0; ix < size_x; ++ix) {
                int64_t i = ix + size_x * (iy + size_y * iz);
                int64_t oi = hr_ratio_x * (ix + size_x * (hr_ratio_y * (iy + size_y * hr_ratio_z * iz)));
                dest[i] = hr_src[oi];
            }
        }
    }
}
