/*------------------------------------------------------------------------------------------------------
  
  File        : AnisotropicTensorDistanceConfidence.h           (GCM Library)

  Description : This class extends the class "AnisotropicTensorDistance" and complete it by adding the statistics of the confidence measure.
  See section 3 and 4.4 of the INRIA Research Report 5845.

  Authors      : Emmanuel Prados (UCLA/INRIA), Christophe Lenglet (INRIA), Jean-Philippe Pons (INRIA)
  
  --------------

  License     : This software is governed by the CeCILL-C license under French law and abiding by the 
  rules of distribution of free software. 

  Users can use, modify and/ or redistribute the software under the terms of the CeCILL-C. In particular, 
  the exercising of this right is conditional upon the obligation to make available to the community the 
  modifications made to the source code of the software so as to contribute to its evolution (e.g. by the 
  mean of the web; i.e. by publishing a web page!).
  
  In this respect, the risks associated with loading, using, modifying and/or developing or reproducing 
  the software by the user are brought to the user's attention, given its Free Software status, which may 
  make it complicated to use, with the result that its use is reserved for developers and experienced 
  professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the 
  suitability of the software as regards their requirements in conditions enabling the security of their
  systems and/or data to be ensured and, more generally, to use and operate it in the same conditions of 
  security. This Agreement may be freely reproduced and published, provided it is not altered, and that 
  no provisions are either added or removed herefrom. 
  
  CeCILL-C FREE SOFTWARE LICENSE AGREEMENT is available in the file 
                          Licence_CeCILL-C_V1-en.txt 
  or at  
                          http://www.cecill.info/index.en.html.

  This Agreement may apply to any or all software for which the holder of the economic rights decides to 
  submit the use thereof to its provisions.
  
  --------------

  Associated publications  : This C++ code corresponds to the implementation of the algorithm presented 
  in the following articles:  
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; 
    "Control Theory and Fast Marching Methods for Brain Connectivity Mapping"; 
    INRIA Research Report 5845 -- UCLA Computer Science Department Technical Report 060004, February 2006.
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; 
    "Control Theory and Fast Marching Methods for Brain Connectivity Mapping";
    Proc. IEEE Computer Society Conference on Computer Vision and Pattern Recognition, New York, NY, I: 1076-1083, June 17-22, 2006.  
  - For more references, we refer to the official web page of the GCM Library and to authors' web pages.

  Please, if you use the GCM library in you work, make sure you will include the reference to the work of 
  the authors in your publications.


  ----------------------------------------------------------------------------------------------------*/


#ifndef ANISOTROPICTENSORDISTANCECONFIDENCE_H
#define ANISOTROPICTENSORDISTANCECONFIDENCE_H

#include "AnisotropicTensorDistance.h"

namespace FastLevelSet {

    template <typename T = float>
    class AnisotropicTensorDistanceConfidence : public AnisotropicTensorDistance<T> {

        protected:
            typedef T Matrix3x3[3][3];

            T *confidenceMean;
            T *confidenceStd;
            T *confidenceMin;

            T *tensor_power;
            T alpha;

        public:
            AnisotropicTensorDistanceConfidence(T *_data, int _width, int _height, int _depth, T* mask, T *tensor, T *_tensor_power, T _alpha, T _dx = T(1), T _dy = T(1), T _dz = T(1)) : AnisotropicTensorDistance<T>(_data,_width,_height,_depth,mask,tensor,_dx,_dy,_dz), tensor_power(_tensor_power), alpha(_alpha) {

                // Allocates the arrays for the statistics of the confidence measure
                confidenceMean = new T[this->size];
                confidenceStd  = new T[this->size];
                confidenceMin  = new T[this->size];

                for (int n=0;n<this->size;n++) {
                    confidenceMean[n] = 0;
                    confidenceStd[n]  = 0;
                    confidenceMin[n]  = FLT_MAX;
                }
            }


            virtual ~AnisotropicTensorDistanceConfidence() {
                delete confidenceMean;
                delete confidenceStd;
                delete confidenceMin;
            }


            T getConfidenceMean(int x, int y, int z) const {
                return confidenceMean[this->_offset(x,y,z)];
            }

            T getConfidenceStd(int x, int y, int z) const {
                return confidenceStd[this->_offset(x,y,z)];
            }

            T getConfidenceMin(int x, int y, int z) const {
                return confidenceMin[this->_offset(x,y,z)];
            }

        protected:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Updates the value at a point
        virtual T _UpdateValue(const int x, const int y, const int z) const {

            // Updates the geodesic distance and the optimal dynamic
            const T val = AnisotropicTensorDistance<T>::_UpdateValue(x,y,z);

            // Gets the optimal dynamic
            const int n = this->_offset(x,y,z);

            const T dyn1 = this->OptDynamics[n][0];
            const T dyn2 = this->OptDynamics[n][1];
            const T dyn3 = this->OptDynamics[n][2];


            T dyn_norm2, dyn_norm;

            if (alpha == 0.0) {

                // Computes its euclidean norm
                dyn_norm2 = dyn1*dyn1 + dyn2*dyn2 + dyn3*dyn3;
                dyn_norm = std::sqrt(dyn_norm2);

            } else if (alpha == -1.0 || alpha > 0) {

                // Computes its alpha norm
                const Matrix3x3 Da = { { tensor_power[n], tensor_power[n+this->size], tensor_power[n+2*this->size] },
                    { tensor_power[n+this->size], tensor_power[n+3*this->size], tensor_power[n+4*this->size] },
                    { tensor_power[n+2*this->size], tensor_power[n+4*this->size], tensor_power[n+5*this->size] }};

                dyn_norm2 = Da[0][0]*dyn1*dyn1 + Da[1][1]*dyn2*dyn2 + Da[2][2]*dyn3*dyn3 + 2*Da[0][1]*dyn1*dyn2 + 2*Da[0][2]*dyn1*dyn3 + 2*Da[1][2]*dyn2*dyn3;
                dyn_norm = std::sqrt(dyn_norm2);

            } else {
                std::cerr << "Unhandled tensor power, exiting..." << std::endl;
                exit(1);
            }

            // Computes the confidence from the neighboring confidence values
            T numMean = dyn_norm;
            T numStd  = dyn_norm2;
            T dyn_norm_y = 0;
            T den = 0;

            if (dyn1>0) {
                den += dyn1;
                numMean += dyn1 * confidenceMean[this->_offset(x-1,y,z)];
                numStd += dyn1 * confidenceStd[this->_offset(x-1,y,z)];
                dyn_norm_y += dyn1 * confidenceMin[this->_offset(x-1,y,z)];
            } else if (dyn1<0) {
                den -= dyn1;
                numMean -= dyn1 * confidenceMean[this->_offset(x+1,y,z)];
                numStd -= dyn1 * confidenceStd[this->_offset(x+1,y,z)];
                dyn_norm_y -= dyn1 * confidenceMin[this->_offset(x+1,y,z)];
            }

            if (dyn2>0) {
                den += dyn2;
                numMean += dyn2 * confidenceMean[this->_offset(x,y-1,z)];
                numStd += dyn2 * confidenceStd[this->_offset(x,y-1,z)];
                dyn_norm_y += dyn2 * confidenceMin[this->_offset(x,y-1,z)];
            } else if (dyn2<0) {
                den -= dyn2;
                numMean -= dyn2 * confidenceMean[this->_offset(x,y+1,z)];
                numStd -= dyn2 * confidenceStd[this->_offset(x,y+1,z)];
                dyn_norm_y -= dyn2 * confidenceMin[this->_offset(x,y+1,z)];
            }

            if (dyn3>0) {
                den += dyn3;
                numMean += dyn3 * confidenceMean[this->_offset(x,y,z-1)];
                numStd += dyn3 * confidenceStd[this->_offset(x,y,z-1)];
                dyn_norm_y += dyn3 * confidenceMin[this->_offset(x,y,z-1)];
            } else if (dyn3<0) {
                den -= dyn3;
                numMean -= dyn3 * confidenceMean[this->_offset(x,y,z+1)];
                numStd -= dyn3 * confidenceStd[this->_offset(x,y,z+1)];
                dyn_norm_y -= dyn3 * confidenceMin[this->_offset(x,y,z+1)];
            }

            if (den != 0) {
                confidenceMean[n] = numMean / den;
                confidenceStd[n] = numStd / den;
                confidenceMin[n] = std::min(dyn_norm, dyn_norm_y / den);
            } else {
//                 std::cerr << "An optimal dynamics with zero magnitude occured!" << std::endl;
                confidenceMean[n] = 0.0;
                confidenceStd[n] = 0.0; // this is a sum of squares, can't be negative!
                confidenceMin[n] = 0.0;
            }

            return val;
        }
    };
}

#endif
