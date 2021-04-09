//   medianfilter.cpp - impelementation of 
//   1D and 2D median filter routines
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#include <memory.h>
#include "Common.hpp"

//   1D MEDIAN FILTER implementation
//     signal - input signal
//     result - output signal
//     NX      - length of the signal
void _medianfilter(const element* signal, element* result, int NX)
{
   //   Move window through all elements of the signal
   for (int i = 2; i < NX - 2; ++i)
   {
      //   Pick up window elements
      element window[5];
      for (int j = 0; j < 5; ++j)
         window[j] = signal[i - 2 + j];
      //   Order elements (only half of them)
      for (int j = 0; j < 3; ++j)
      {
         //   Find position of minimum element
         int min = j;
         for (int k = j + 1; k < 5; ++k)
            if (window[k] < window[min])
               min = k;
         //   Put found minimum element in its place
         const element temp = window[j];
         window[j] = window[min];
         window[min] = temp;
      }
      //   Get result - the middle element
      result[i - 2] = window[2];
   }
}

//   1D MEDIAN FILTER wrapper
//     signal - input signal
//     result - output signal
//     NX      - length of the signal
void medianfilter(element* signal, element* result, int NX)
{
   //   Check arguments
   if (!signal || NX < 1)
      return;
   //   Treat special case NX = 1
   if (NX == 1)
   {
      if (result)
         result[0] = signal[0];
      return;
   }
   //   Allocate memory for signal extension
   element* extension = new element[NX + 4];
   //   Check memory allocation
   if (!extension)
      return;
   //   Create signal extension
   memcpy(extension + 2, signal, NX * sizeof(element));
   for (int i = 0; i < 2; ++i)
   {
      extension[i] = signal[1 - i];
      extension[NX + 2 + i] = signal[NX - 1 - i];
   }
   //   Call median filter implementation
   _medianfilter(extension, result ? result : signal, NX + 4);
   //   Free memory
   delete[] extension;
}

//   2D MEDIAN FILTER implementation
//     image  - input image
//     result - output image
//     NX      - width of the image
//     NY      - height of the image
void _medianfilter(const element* image, element* result, int NX, int NY, int ksize)
{
   int ksize2 = ksize * ksize;
   int kBoundary = (ksize - 1) / 2;
   //   Move window through all elements of the image
   for (int iY = kBoundary; iY < NY - kBoundary; ++iY)
      for (int iX = kBoundary; iX < NX - kBoundary; ++iX)
      {
         //   Pick up window elements
         int k = 0;
         element window[ksize2];
         for (int j = iY - kBoundary; j < iY + kBoundary + 1; ++j)
            for (int i = iX - kBoundary; i < iX + kBoundary + 1; ++i)
               window[k++] = image[j * NX + i];
         //   Order elements (only half of them)
         for (int j = 0; j < (ksize2+1)/2; ++j)
         {
            //   Find position of minimum element
            int min = j;
            for (int l = j + 1; l < ksize2; ++l)
            if (window[l] < window[min])
               min = l;
            //   Put found minimum element in its place
            const element temp = window[j];
            window[j] = window[min];
            window[min] = temp;
         }
         //   Get result - the middle element
         result[(iY - kBoundary) * (NX - 2*kBoundary) + iX - kBoundary] = window[(ksize2-1)/2];
      }
}

//   2D MEDIAN FILTER wrapper
//     image  - input image
//     result - output image
//     NX      - width of the image
//     NY      - height of the image
//     ksize   - window size = ksize*ksize
void medianfilter(element* image, element* result, int NX, int NY, int ksize, int fillType)
{    
   int kBoundary = (ksize - 1) / 2;
   //   Check arguments
   if (!image || NX < 1 || NY < 1)
      return;
   //   Allocate memory for signal extension
   element* extension = new element[(NX + 2*kBoundary) * (NY + 2*kBoundary)];
   //   Check memory allocation
   if (!extension)
      return;
   //   Create image extension
   for (int i = 0; i < NY; ++i)
   {
      memcpy(extension + (NX + 2*kBoundary) * (i + kBoundary) + kBoundary, image + NX * i, NX * sizeof(element));
      if (fillType == 1) {
        for (int j = 0; j < kBoundary; ++j) {
            extension[(NX + 2*kBoundary) * (i + kBoundary) + j] = image[NX * i];
            extension[(NX + 2*kBoundary) * (i + kBoundary + 1) - j -1] = image[NX * (i + 1) - 1];
        }
      } else {
            for (int j = 0; j < kBoundary; ++j) {
                extension[(NX + 2*kBoundary) * (i + kBoundary) + j] = 0;
                extension[(NX + 2*kBoundary) * (i + kBoundary + 1) - j -1] = 0;
            }
      }
   }
   for (int j = 0; j < kBoundary; ++j)
   {
       if (fillType == 1) {
            //   Fill first line of image extension
            memcpy(extension + (NX + 2*kBoundary) * j, extension + (NX + 2*kBoundary) * kBoundary, (NX + 2*kBoundary) * sizeof(element));
            //   Fill last line of image extension
            memcpy(extension + (NX + 2*kBoundary) * (NY + kBoundary + j), extension + (NX + 2*kBoundary) * (NY + kBoundary - 1), (NX + 2*kBoundary) * sizeof(element));
      } else {
            for (int i = 0; i < NX + 2*kBoundary; ++i) {
                //   Fill first line of image extension
                extension[(NX + 2*kBoundary) * j + i] = 0;
                //   Fill last line of image extension
                extension[(NX + 2*kBoundary) * (NY + kBoundary + j) + i] = 0;
            }
       }
   }
   //   Call median filter implementation
   _medianfilter(extension, result ? result : image, NX + 2*kBoundary, NY + 2*kBoundary, ksize);
   //   Free memory
   delete[] extension;
}

/*! \brief Apply a median filter to filter the extrame value of 2D gradient vector
 */
// template<typename ValueType>
// void applyMedianFilterTo2DVector(scai::lama::DenseVector<ValueType> &vecter2D, scai::IndexType NX, scai::IndexType NY, scai::IndexType spatialFDorder)
// {    
//     element signal_2D[ NY * NX ] = {0};
//     for (int i = 0; i < NY * NX; i++) {
//         signal_2D[i] = vecter2D.getValue(i);
//     }
//     element *result;
//     result = (element *)malloc(NY * NX * sizeof(element));
// 
//     medianfilter(signal_2D, result, NX, NY, spatialFDorder+3, 0);
//     
//     for (int i = 0; i < NY * NX; i++) {
//         vecter2D.setValue(i, *(result + i));
//     }
// }
