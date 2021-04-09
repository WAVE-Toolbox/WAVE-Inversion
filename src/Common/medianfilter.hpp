//   medianfilter.hpp - declarations for 
//   1D and 2D median filter routines
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

// #ifndef _MEDIANFILTER_H_
// #define _MEDIANFILTER_H_
// 
// //   Signal/image element type
// typedef float element;
// 
// //   1D MEDIAN FILTER, window size 5
// //     signal - input signal
// //     result - output signal, NULL for inplace processing
// //     NX      - length of the signal
// void medianfilter(element* signal, element* result, int NX) = 0;
// 
// //   2D MEDIAN FILTER, window size 3x3
// //     image  - input image
// //     result - output image, NULL for inplace processing
// //     NX      - width of the image
// //     NY      - height of the image
// void medianfilter(element* image, element* result, int NX, int NY, int ksize, int fillType) = 0;
// 
// #endif
