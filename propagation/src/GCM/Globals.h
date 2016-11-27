/*------------------------------------------------------------------------------------------------------
  
  File        : Globals.h   (GCM Library)

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


#ifndef FASTLEVELSETGLOBALS_H
#define FASTLEVELSETGLOBALS_H

#include <list>

#ifndef sqr
#define sqr(x) ((x)*(x))
#endif

namespace FastLevelSet {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constants
    const float DELTA = 1.0f;
    const float EPS = 1e-6f;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Types
    struct stCoord {
        public:
        int x,y,z,n;
        stCoord(int _x, int _y, int _z, int _n) : x(_x), y(_y), z(_z), n(_n) {}
    };

    typedef std::list<stCoord> CoordList;
    typedef CoordList::iterator CoordListIterator;
    typedef CoordList::const_iterator const_CoordListIterator;
}

#endif
