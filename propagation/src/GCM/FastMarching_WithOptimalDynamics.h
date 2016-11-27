/*------------------------------------------------------------------------------------------------------
  
  File        : FastMarching_WithOptimalDynamics.h   (GCM Library)

  Description : This class complete the "FastMarching" class by adding the necessary tools for dealing 
                with the Optimal Dynamics.

  Authors      : Emmanuel Prados (UCLA/INRIA), Christophe Lenglet (INRIA), Jean-Philippe Pons (INRIA)
  
  --------------

  License     : This software is governed by the CeCILL-C license under French law and abiding by the 
  rules of distribution of free software. 

  Users can use, modify and/ or redistribute the software under the terms of the CeCILL-C. In particular, 
  the exercising of this right is conditional upon the obligation to make available to the community the 
  modifications made to the source code of the software so as to contribute to its evolution (e.g. by the 
  mean of the web; i.e. by publishing a web page).
  
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
    Control Theory and Fast Marching Methods for Brain Connectivity Mapping; INRIA Research Report 5845 
    UCLA Computer Science Department Technical Report 060004, February 2006.
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; 
    Control Theory and Fast Marching Methods for Brain Connectivity Mapping; 
    Proc. IEEE Computer Society Conference on Computer Vision and Pattern Recognition, New York, NY, 
    I: 1076-1083, June 17-22, 2006.  
  - For more references, we refer to the official web page of the GCM Library and to authors' web pages.

  Please, if you use the GCM library in you work, make sure you will include the reference to the work 
  of the authors in your publications.
  
  
----------------------------------------------------------------------------------------------------*/



#ifndef FASTLEVELSETFASTMARCHING_WITHOPTIMALDYNAMICS_H
#define FASTLEVELSETFASTMARCHING_WITHOPTIMALDYNAMICS_H

#include "FastMarching.h"

namespace FastLevelSet {

/*======================================================================================================

     The optimal Dymamics are not usefull if we only want to compute the viscosity solution of the 
     considered equation by the Fast Marching Method. Nevertheless optimal Dymamics is usefull in
     many applications as for example fibers tracking in DTI (see papers associated to this code).

======================================================================================================*/


    template <typename T = float>
    class FastMarching_WithOptDynamics : public FastMarching<T,1> {

    typedef T Vector3[3];          // 3D-vector

    protected:

        Vector3 *OptDynamics;       // Array containing the Optimal Dynamics (array of Vector3):
                                    // Theoretically the dynamics is a vector whose
                                    // the dimension depends on the dimension of the space variable.
                                    // if width!=0, height==0 and depth==0 then dim of dynamics is one.
                                    // if width!=0, height!=0 and depth==0 then dim of dynamics is two.
                                    // In practice, we use a 3D-vector in order to deal indifferently
                                    // with any dimension <=3...

    public:
        // Constructor
        FastMarching_WithOptDynamics(T *_data, int _width, int _height, int _depth) : FastMarching<T,1>(_data,_width,_height,_depth) {

            // Allocation and initialization of the array "OptDynamics":
            OptDynamics = new Vector3[this->size];

            for (int n=0;n<this->size;n++) {
                OptDynamics[n][0] = 0;
                OptDynamics[n][1] = 0;
                OptDynamics[n][2] = 0;
            }
        }

        // Destructor
        virtual ~FastMarching_WithOptDynamics() {
            delete [] OptDynamics;
        }

        // accessors & Co
        T getOptDynamics(int x, int y, int z, int component) const {
            const int n = this->_offset(x,y,z);
            return OptDynamics[n][component];
        }

        virtual void setOptDynamics(int x, int y, int z, T component1, T component2, T component3) const {
            const int n = this->_offset(x,y,z);
            OptDynamics[n][0]= component1;
            OptDynamics[n][1]= component2;
            OptDynamics[n][2]= component3;
        }
    };  // end of class "FastMarching_WithOptDynamics"
}  // end of namespace "FastLevelSet"

#endif
