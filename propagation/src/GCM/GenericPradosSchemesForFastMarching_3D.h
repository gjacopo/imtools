/*------------------------------------------------------------------------------------------------------
  
  File        : GenericPradosSchemesForFastMarching_3D.h   (GCM Library)

  Authors      : Emmanuel Prados (UCLA/INRIA), Christophe Lenglet (INRIA), Jean-Philippe Pons (INRIA)

  Description : This method explicits the scheme proposed by Prados etal. in dimension 3 (this scheme is
                shortly described in sections 4.3 of INRIA Research Report 5845  -- in particular, see 
		sections 4.3.1 and 4.3.2).
  
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


#ifndef GENERICPRADOSSCHEMESFORFASTMARCHING_3D_H
#define GENERICPRADOSSCHEMESFORFASTMARCHING_3D_H

#include "Globals.h"
#include "FastMarching_WithOptimalDynamics.h"

/**************************************************/
// Dimension 3 =====================================
/**************************************************/
/*

Here, we deal with the equation of the form:

$$ H(x,\nabla u(x)) = 0 . $$

Let $H^*$ be the Legendre transform of $H$.
For $i \in [1..N]$, we denote

$$ H_i(x,p) = sup_{a\in Dom(H^*), a_i = 0}  a \cdot p - H^*(x,a)$$.
$$ H_ij(x,p) = sup_{a\in Dom(H^*), a_i = 0 and a_j = 0,  j \neq i}  a \cdot p - H^*(x,a)$$.

*/
/**************************************************/


// Note about the optimal dynamics: 
// The preservation and the transmission of the optimal dymamics is not usefull if we only want 
// to compute the viscosity solution of the considered equation by the Fast Marching Method. 
// Nevertheless the computation of the optimal dymamics (function $f$ in Prados's papers) is 
// necessary inside (i.e. when we develop) the functions
//      bool eqSolverOnPart_with_s1s2s3_nonNull(...)
//      bool eqSolverOnPart_withOne_si_Null(...)
//      bool eqSolverOnPart_withTwo_si_Null(...)
// in order to be able to know if the considered simplex is a good candidate.
// Note: preservation and transmission of the optimal dymamics can be usefull in many applications
// as for example fibers tracking in DTI.


namespace FastLevelSet {

    template <typename T = float>
    class PradosSchemesForFastMarching_3D : public FastMarching_WithOptDynamics<T> {

    public:
        // Constructor
        PradosSchemesForFastMarching_3D(T *_data, int _width, int _height, int _depth) : FastMarching_WithOptDynamics<T>(_data,_width,_height,_depth) {}

        // Destructor
        virtual ~PradosSchemesForFastMarching_3D() {}

    protected:

        //////////////////////////////////////////////////////
        // These next THREE Functions must be overcharged !
        //////////////////////////////////////////////////////


        //////////////////////////////////////////////////////////////////////////
        //
        //  Let us denote:
        //
        //      $$ H_{s1s2s3}(x,p) = sup_{a\in Dom_{s1s2s3}}  a \cdot p - H^*(x,a)$$.
        //      where
        //      $$ Dom_{s1s2s3}  = \{ a\in Dom(H^*) such that for all i, sign(a_i)*si < 0,
        //
        //  eqSolverOnPart_with_s1s2s3_nonNull must solve the following equation (in $t$) :
        //
        //  $ H_{s1s2s3}(x,p_t) = 0 $     (equation (1)),
        //
        //  where $[p_t]_i = [ t-u(x+sihiei)] / (-sihi)$.
        //
        //  Solve equation (1) is equivalent to solve equation
        //
        //  $ H(x,p_t) = 0 $     (equation (2)),
        //
        //  and then, amongst the solutions of equation (2),  chose
        //  the solution $t_0$ of such that:
        //  $\nabla H (x,p_{t_0}) * si <0 $.
        //
        //  Notes:
        //  1)  $\nabla H (x,p_{t_0})$ corresponds with the optimal
        //      control of $ H(x,p_{t_0}) = 0 $.
        //
        //  2)  Here, the $hi$ are not given as parameters, since they
        //  are supposed constant, and known in the next inherited
        //  classes...
        //
        //  3) the parameter Ui corresponds with the value $u(x+sihiei)$
        //
        //////////////////////////////////////////////////////////////////////////

        virtual bool eqSolverOnPart_with_s1s2s3_nonNull(
            const T U1, const T U2, const T U3,                 // Values of the solution at the considered neigborhood voxels,
            const int s1, const int s2, const int s3,           // signs associated to the considered sector,
            const int x, const int y, const int z,              // coordinates of the considered voxel,
            T &Root,                                            // solution.
            T &optDymamics1, T &optDymamics2, T &optDymamics3   // optimal dynamic associated to the solution.
            ) const = 0;

        virtual bool eqSolverOnPart_withOne_si_Null(
            const T U1, const T U2, const T U3,
            const int s1, const int s2, const int s3,
            const int x, const int y, const int z,
            const int indice_si_EqualZero,
            T &Root,
            T &optDymamics1,    T &optDymamics2,    T &optDymamics3
            ) const = 0;

        virtual bool eqSolverOnPart_withTwo_si_Null(
            const T U1, const T U2, const T U3,
            const int s1, const int s2, const int s3,
            const int x, const int y, const int z,
            const int indice_si_DiffZero,
            T &Root,
            T &optDymamics1,    T &optDymamics2,    T &optDymamics3
            ) const = 0;


        //////////////////////////////////////////////////////
        // The next function must NOT be overcharged !
        //////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Update of a point
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////

        virtual T _UpdateValue(const int x, const int y, const int z) const {

            // Let us deal with the boundary:
            if ((x==0) || (x==this->width-1)  ||
                (y==0) || (y==this->height-1) ||
                (z==0) || (z==this->depth-1))
                return this->big;

            // Initialisation of s1,s2,S3.
            int s1=0, s2=0, s3=0;

            // Initialisation of DoExistSol.
            bool DoExistSol[3][3][3];
            for(s1=0; s1<3; s1++) for(s2=0; s2<3; s2++) for(s3=0; s3<3; s3++) {
                DoExistSol[s1][s2][s3]=false;
            }

            // Initialisation of Root.
            T Root[3][3][3];
            for(s1=0; s1<3; s1++) for(s2=0; s2<3; s2++) for(s3=0; s3<3; s3++) {
                Root[s1][s2][s3] = this->big;
            }

            // Initialisation of the optimal dynamics.
            // In dimension 3, the dynamics (function f associated to an HJB eq.) is a 3D-vector.
            // optDymamics1, optDymamics2, optDymamics3 represente the its componantes.
            // As for "Root" and "DoExistSol", we define and compute an optimal
            // Dymamics on each simplex and each "boundaries" of simplexes.

            T optDymamics1[3][3][3];
            T optDymamics2[3][3][3];
            T optDymamics3[3][3][3];

            for(s1=0; s1<3; s1++) for(s2=0; s2<3; s2++) for(s3=0; s3<3; s3++) {
                optDymamics1[s1][s2][s3] = 0;
                optDymamics2[s1][s2][s3] = 0;
                optDymamics3[s1][s2][s3] = 0;
            }


            // ===========================================

            // Computation of Roots such that opt control
            // is in interior of the various Ds1s2s3:

            for(s1=-1; s1<=1; s1+=2)
                for(s2=-1; s2<=1; s2+=2)
                    for(s3=-1; s3<=1; s3+=2) {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_with_s1s2s3_nonNull(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }

            // ===========================================

            // Computation of Roots such that opt control
            // is in interior of the intersection of Ds1s2s3
            // and set of a such that f_i(x,a) = 0

            s1=0;
            for(s2=-1; s2<=1; s2+=2)
                for(s3=-1; s3<=1; s3+=2) {
                    if ( ! (DoExistSol[-1+1][s2+1][s3+1] || DoExistSol[1+1][s2+1][s3+1]) ) {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withOne_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                1,  // Indice_i_Such_siEqualTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
                    else{
                        DoExistSol[s1+1][s2+1][s3+1] = true;  // useful for the next step!
                    }
                }

            s2=0;
            for(s1=-1; s1<=1; s1+=2)
                for(s3=-1; s3<=1; s3+=2) {
                    if ( ! (DoExistSol[s1+1][-1+1][s3+1] || DoExistSol[s1+1][1+1][s3+1]) ) {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withOne_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                2,  // Indice_i_Such_siEqualTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
                    else{
                        DoExistSol[s1+1][s2+1][s3+1] = true;  // useful for the next step!
                    }
                }

            s3=0;
            for(s2=-1; s2<=1; s2+=2)
                for(s1=-1; s1<=1; s1+=2) {
                    if ( ! (DoExistSol[s1+1][s2+1][-1+1] || DoExistSol[s1+1][s2+1][1+1]) ) {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withOne_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                3,  // Indice_i_Such_siEqualTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
                    else {
                        DoExistSol[s1+1][s2+1][s3+1] = true;  // useful for the next step!
                    }
                }


            // ===========================================
            // computation of Roots such that opt control
            // is in Interior of the intersection of Ds1s2s3
            // and set of a such that f_i(x,a) = 0

            s1=0;
            s2=0;
            for(s3=-1; s3<=1; s3+=2) {
                // More comments could be required here! ....
                if ( ! (
                    DoExistSol[ 0+1][-1+1][s3+1]
                    || DoExistSol[ 0+1][ 1+1][s3+1]
                    || DoExistSol[-1+1][ 0+1][s3+1]
                    || DoExistSol[ 1+1][ 0+1][s3+1]
                    )
                    )
                    {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withTwo_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                3,  // Indice_i_Such_siDifferentTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
            }

            s1=0;
            s3=0;
            for(s2=-1; s2<=1; s2+=2) {
                // More comments could be required here! ....
                if ( ! (
                    DoExistSol[ 0+1][s2+1][-1+1]
                    || DoExistSol[ 0+1][s2+1][ 1+1]
                    || DoExistSol[-1+1][s2+1][ 0+1]
                    || DoExistSol[ 1+1][s2+1][ 0+1]
                    )
                    )
                    {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withTwo_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                2,  // Indice_i_Such_siDifferentTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
            }

            s2=0;
            s3=0;
            for(s1=-1; s1<=1; s1+=2) {
                // More comments could be required here! ....
                if ( ! (
                    DoExistSol[s1+1][ 0+1][-1+1]
                    || DoExistSol[s1+1][ 0+1][ 1+1]
                    || DoExistSol[s1+1][-1+1][ 0+1]
                    || DoExistSol[s1+1][ 1+1][ 0+1]
                    )
                    )
                    {
                        DoExistSol[s1+1][s2+1][s3+1]
                            = eqSolverOnPart_withTwo_si_Null(
                                this->_GetValue(x+s1,y,z),
                                this->_GetValue(x,y+s2,z),
                                this->_GetValue(x,y,z+s3),
                                s1,s2,s3,
                                x,y,z,
                                1,  // Indice_i_Such_siDifferentTo0;
                                Root[s1+1][s2+1][s3+1],
                                optDymamics1[s1+1][s2+1][s3+1],optDymamics2[s1+1][s2+1][s3+1], optDymamics3[s1+1][s2+1][s3+1]
                                );
                    }
            }


            // ===========================================================================


            // The solution is then the smallest Root:
            T minRoot = this->big;      //  the smallest root
            int gs1=0, gs2=0, gs3=0;    //  give the simplex of the smallest
                                        //  root (usefull for recovering the optimal dynamics).

            for(s1=0; s1<3; s1++)
                for(s2=0; s2<3; s2++)
                    for(s3=0; s3<3; s3++)
                        if (minRoot > Root[s1][s2][s3])
                        {
                            minRoot = Root[s1][s2][s3];     // the current min is then Root[s1][s2][s3]
                            gs1=s1; gs2=s2; gs3=s3;         // remember the simplex associated to the current min
                        };


            // Let us record the optimal Dymamics associated to the adequate simplex.
            // The packed optimal dynamics is the one associated to the used simplex.
            // Note: Preservation of the optimal Dymamics is not usefull if we only
            // want to compute the viscosity solution of the considered equation by
            // the Fast Marching Method. Nevertheless optimal Dymamics is usefull in
            // many application as for example fibers tracking in DTI (see papers
            // associated to this code).

            setOptDynamics(x,y,z, optDymamics1[gs1][gs2][gs3],optDymamics2[gs1][gs2][gs3],optDymamics3[gs1][gs2][gs3]);

            // we can then return the value for update...
            return minRoot;

        }// End of the methode "_UpdateValue".

    };  // End of the definition of the class.

} // End of namespace FastLevelSet.

#endif
