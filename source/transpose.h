#pragma once


template <typename T>
void transpose(size_t width1, size_t width2,
    const T* pSrc, int srcLineStride,
    T* pDst, int dstLineStride) noexcept
{
    if ((width1 == 0) || (width2 == 0) || (pSrc == nullptr) || (pDst == nullptr) || (srcLineStride < width1) || (dstLineStride < width2)) {
        return;
    }


    const size_t bw = 256;
    const size_t bh = 64;
    const size_t bx_end = (width1 / bw) * bw;
    const size_t by_end = (width2 / bh) * bh;
    for (size_t by = 0; by < by_end; by += bh) {
        for (size_t bx = 0; bx < bx_end; bx += bw) {
            auto offset_s = by * srcLineStride + bx;
            auto offset_d = bx * dstLineStride + by;

#pragma _NEC ivdep
#pragma ivdep
            for (size_t iy = 0; iy < bh; ++iy) {
#pragma vector
                for (size_t ix = 0; ix < bw; ++ix) {
                    auto src = pSrc[offset_s + iy * srcLineStride + ix];
                    pDst[offset_d + ix * dstLineStride + iy] = src;
                }
            }
        }
    }

    const size_t widthRemain = width1 - bx_end;
    const size_t heightRemain = width2 - by_end;
    {
        size_t by = by_end;
        for (size_t bx = 0; bx < bx_end; bx += bw) {
            auto offset_s = by * srcLineStride + bx;
            auto offset_d = bx * dstLineStride + by;

#pragma _NEC ivdep
#pragma ivdep
            for (size_t iy = 0; iy < heightRemain; ++iy) {
#pragma vector
                for (size_t ix = 0; ix < bw; ++ix) {
                    auto src = pSrc[offset_s + iy * srcLineStride + ix];
                    pDst[offset_d + ix * dstLineStride + iy] = src;
                }
            }
        }
    }
    {
        for (size_t by = 0; by < by_end; by += bh) {
            size_t bx = bx_end;
            auto offset_s = by * srcLineStride + bx;
            auto offset_d = bx * dstLineStride + by;

#pragma _NEC ivdep
#pragma ivdep
            for (size_t iy = 0; iy < bh; ++iy) {
#pragma vector
                for (size_t ix = 0; ix < widthRemain; ++ix) {
                    auto src = pSrc[offset_s + iy * srcLineStride + ix];
                    pDst[offset_d + ix * dstLineStride + iy] = src;
                }
            }
        }
    }
    {
        size_t by = by_end;
        size_t bx = bx_end;
        auto offset_s = by * srcLineStride + bx;
        auto offset_d = bx * dstLineStride + by;

#pragma _NEC ivdep
#pragma ivdep
        for (size_t iy = 0; iy < heightRemain; ++iy) {
            for (size_t ix = 0; ix < widthRemain; ++ix) {
                auto src = pSrc[offset_s + iy * srcLineStride + ix];
                pDst[offset_d + ix * dstLineStride + iy] = src;
            }
        }
    }
}
