/*------------------------------------------------------------------------------------------------------
  
  File        : AnisotropicTensorDistance.h     (GCM Library)

  Description : This class specifies the virtual methods of the class "GenericPradosSchemesForFastMarching_3D" for the special case of the 3D anisotropic Eikonal equation.

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

  --------------
  
  Technical Comments and Detailled description:

  - This class specifies the virtual methods of the class "GenericPradosSchemesForFastMarching_3D" 
  for the special case of the 3D anisotropic Eikonal equation.
  The functions we overcharge are:
   * virtual bool eqSolverOnPart_with_s1s2s3_nonNul(...)
   * virtual bool eqSolverOnPart_withOne_si_Nul(...)
   * virtual bool eqSolverOnPart_withTwo_si_Nul(...)

  $(s_1,s_2,s_3)$ being fixed,  these methods allow to solve the equations   
  $$ sup_{a in B} g_s(a,t) = 0 $$
  where $g_s(a,t) = - f(x; a).P_{x,s,U}(t)-l(x; a)$ 
  and where, according to the different cases, 
  * $B$ is the interior set of the $As$   (eqSolverOnPart_with_s1s2s3_nonNu)
  * $B$ is the union of the 2D facets of $As$  (eqSolverOnPart_withOne_si_Nul)
  * $B$ is the union of the 1D edges of $As$  (eqSolverOnPart_withTwo_si_Nul)
  See section 4.3.2 of INRIA RR-5845 (page 11).
  Let us note that these methods compute the $t_s$ and the associated dynamics.

  - The goal of the class being to specify the virtual methods of the class  
  "GenericPradosSchemesForFastMarching_3D" for the special case of the 3D anisotropic Eikonal equation,
  we then need to describe the data defining the Hamiltonian.
  In this class, we then need to declare and construct all the associated structures.

  - Note about the optimal dynamics:
  The optimal dymamics (so its preservation and its transmission in the paramters) ARE NOT USEFULL IF we
  only want to compute the viscosity solution of the considered equation by the Fast Marching Method. 
  Nevertheless the computation of the optimal dymamics is necessary in order to be able to chose the good 
  simplex. The Preservation and the transmission of the optimal Dymamics can be usefull in many applications 
  as for example fibers tracking in DTI (see INRIA RR-5845).


  ----------------------------------------------------------------------------------------------------*/




#ifndef ANISOTROPICTENSORDISTANCE_H
#define ANISOTROPICTENSORDISTANCE_H

#include <iostream>
#include <cmath>
#include <cassert>
#include "GenericPradosSchemesForFastMarching_3D.h"

/**************************************************/
//
//  This class (AnisotropicTensorDistance) inherits of
//  the class: "PradosSchemesForFastMarching_3D"
//
//  Here we have to overcharge the functions
//   - virtual bool eqSolverOnPart_with_s1s2s3_nonNul(...)
//   - virtual bool eqSolverOnPart_withOne_si_Nul(...)
//   - virtual bool eqSolverOnPart_withTwo_si_Nul(...)
//
/**************************************************/




namespace FastLevelSet {


/**************************************************
  The class AnisotropicTensorDistance inherits of
  the class: "PradosSchemesForFastMarching_3D"
**************************************************/

    template <typename T = float>
    class AnisotropicTensorDistance : public PradosSchemesForFastMarching_3D<T> {

        /***************************************************************
           Let us declare all the structures describing the data
           charaterizing the 3D Eikonal equation considered.
           We also declare some variables required to describe the scheme.
        ****************************************************************/


    protected:
    T dx,dy,dz;                               // voxel sizes
        typedef T Vector3[3];
        typedef T Matrix3x3[3][3];
        typedef T Matrix2x2[2][2];

        // we start with the declaration of the class "Hamiltonian".

          /***************************************************************
            The Hamiltonian that we have implemented here is
                           $  |A_x p |^2 - 1  $.
            The Hamiltonian is equivalent to the Hamiltonian
                     $  H_{AEik}(x,p) = |A_x p | - 1 $
            given at page 7 of the  INRIA report RR-5845.
            The $ H_i(x,p) $  and $ H_{ij}(x,p) $ described at page 12 (section 4.3.3) of the INRIA report RR-5845
            are, in fact, the ones associated to $|A_x p |^2 - 1$.
            There are no differences, except that we need to renormalize the optimal dynamics. 
            We will explain better this point in our forthcomming journal paper...
          ****************************************************************/

        class Hamiltonian {
        public:
            Matrix3x3 B, invB;
            Matrix2x2 invinvB1, invinvB2, invinvB3;
            T invsqinvA1, invsqinvA2, invsqinvA3;

            // Initialization from B
            bool InitFromB(const Matrix3x3 &_B) {

                // Storage of B
                CopyMatrix3x3(_B,B);

                // Computation of the inverse of B
                // if (!InverseSymmetricMatrix3x3(B,invB)) return false;
                if (Matrix3x3IsZero(B)) return false;
                else PseudoInverseSymmetricMatrix3x3(B,invB);

                // Computation of invinvB_i: the inverse of invB after removing the i^th row and i^th column
                for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
                    invinvB1[i][j] = invB[i+1][j+1];
                    invinvB2[i][j] = invB[2*i][2*j];
                    invinvB3[i][j] = invB[i][j];
                }

                if (Matrix2x2IsZero(invinvB1)) return false;
                else PseudoInverseMatrix2x2(invinvB1);

                if (Matrix2x2IsZero(invinvB2)) return false;
                else PseudoInverseMatrix2x2(invinvB2);

                if (Matrix2x2IsZero(invinvB3)) return false;
                else PseudoInverseMatrix2x2(invinvB3);

//                 if (!InverseMatrix2x2(invinvB1) || !InverseMatrix2x2(invinvB2) || !InverseMatrix2x2(invinvB3)) return false;

                // Computation of invinvAi: the inverse of the square norm of the i^th row of invA
                invsqinvA1 = 1/invB[0][0];
                invsqinvA2 = 1/invB[1][1];
                invsqinvA3 = 1/invB[2][2];

                return true;
            }

            // Initialization from A
            bool InitFromA(const Matrix3x3 &_A) {
                // TODO
                return false;
            }

        private:
            // Matrix operations
            void DisplayMatrix3x3(const Matrix3x3 &M) {
                for (int i=0;i<3;i++)
                    std::cout << M[i][0] << " " << M[i][1] << " " << M[i][2] << std::endl;
                std::cout << std::endl;
            }

            void DisplayMatrix2x2(const Matrix2x2 &M) {
                for (int i=0;i<2;i++)
                    std::cout << M[i][0] << " " << M[i][1]  << std::endl;
                std::cout << std::endl;
            }

            void MultiplyMatrix3x3(const Matrix3x3 &A, const Matrix3x3 &B, Matrix3x3 &R) {
                for (int i=0; i<3; ++i) {
                    for (int j=0; j<3; ++j) {
                        R[i][j] = 0.0;
                        for (int k=0; k<3; ++k) R[i][j] += A[i][k] * B[k][j];
                    }
                }
            }

            void MultiplyMatrix2x2(const Matrix2x2 &A, const Matrix2x2 &B, Matrix2x2 &R) {
                for (int i=0; i<2; ++i) {
                    for (int j=0; j<2; ++j) {
                        R[i][j] = 0.0;
                        for (int k=0; k<2; ++k) R[i][j] += A[i][k] * B[k][j];
                    }
                }
            }

            void CopyMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {
                for (int i=0;i<3;i++) for (int j=0;j<3;j++) out[i][j] = in[i][j];
            }

            void CopyMatrix2x2(const Matrix2x2 &in, Matrix2x2 &out) {
                for (int i=0;i<2;i++) for (int j=0;j<2;j++) out[i][j] = in[i][j];
            }

            void TransposeMatrix2x2(const Matrix2x2 &in, Matrix2x2 &out) {
                for (int i=0;i<2;i++) for (int j=0;j<2;j++) out[i][j] = in[j][i];
            }

            void TransposeMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {
                for (int i=0;i<3;i++) for (int j=0;j<3;j++) out[i][j] = in[j][i];
            }

            bool Matrix2x2IsZero(const Matrix2x2 &M) {
                return (M[0][0]==0 && M[0][1]==0 && M[1][1]==0);
            }

            bool Matrix3x3IsZero(const Matrix3x3 &M) {
                return (M[0][0]==0 && M[0][1]==0 && M[0][2]==0 && M[1][1]==0 && M[1][2]==0 && M[2][2]==0);
            }

            bool InverseUpperTriangularMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {
                if (in[0][0] == 0 || in[1][1] == 0 || in[2][2] == 0) return false;
                out[0][0] = 1/in[0][0];
                out[0][1] = -in[0][1]/in[0][0]/in[1][1];
                out[0][2] = (in[0][1]*in[1][2] - in[0][2]*in[1][1])/in[0][0]/in[1][1]/in[2][2];
                out[1][0] = 0;
                out[1][1] = 1/in[1][1];
                out[1][2] = -in[1][2]/in[1][1]/in[2][2];
                out[2][0] = 0;
                out[2][1] = 0;
                out[2][2] = 1/in[2][2];
                return true;
            }

            bool InverseSymmetricMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {

                // for computing the determinant,
                // we develop with respect to the first row:
                T firstCrossProduct  = in[1][1]*in[2][2]-in[1][2]*in[1][2];
                T secondCrossProduct = in[0][1]*in[2][2]-in[0][2]*in[1][2];
                T thirdCrossProduct  = in[0][1]*in[1][2]-in[0][2]*in[1][1];

                T determinant =
                        in[0][0]*firstCrossProduct
                    -   in[0][1]*secondCrossProduct
                    +   in[0][2]*thirdCrossProduct;

                if (determinant == 0) return false;

                out[0][0] =  (  firstCrossProduct   )/determinant;
                out[0][1] =  ( - secondCrossProduct )/determinant;
                out[0][2] =  (  thirdCrossProduct   )/determinant;

                out[1][0] = out[0][1];
                out[1][1] =  (-in[0][2]*in[0][2]+in[0][0]*in[2][2])/determinant;
                out[1][2] =  (-in[0][0]*in[1][2]+in[0][2]*in[0][1])/determinant;

                out[2][0] = out[0][2];
                out[2][1] = out[1][2];
                out[2][2] =  (-in[0][1]*in[0][1]+in[0][0]*in[1][1])/determinant;

                return true;
            }

            void PseudoInverseSymmetricMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {
                Matrix3x3 M,N;
                MultiplyMatrix3x3(in,in,M);
                InverseSymmetricMatrix3x3(M,N);
                MultiplyMatrix3x3(N,in,out);
            }

            //bool CholeskyMatrix3x3(const Matrix3x3 &in, Matrix3x3 &out) {
            //  if (in[0][0]<0) return false;
            //  out[0][0] = (T)std::sqrt(in[0][0]);
            //  out[0][1] = in[0][1] / out[0][0];
            //  out[0][2] = in[0][2] / out[0][0];
            //  out[1][0] = 0;
            //  const T d2 = in[1][1] - out[0][1]*out[0][1];
            //  if (d2<0) return false;
            //  out[1][1] = (T)std::sqrt(d2);
            //  out[1][2] = ( in[1][2] - out[0][1]*out[0][2] ) / out[1][1];
            //  out[2][0] = 0;
            //  out[2][1] = 0;
            //  const T d3 = in[2][2] - out[0][2]*out[0][2] - out[1][2]*out[1][2];
            //  if (d3<0) return false;
            //  out[2][2] = (T)std::sqrt(d3);
            //  return true;
            //}

            bool InverseMatrix2x2(Matrix2x2 &M) {
                const T det = M[0][0]*M[1][1]-M[0][1]*M[1][0];
                if (det==0) return false;
                std::swap(M[0][0],M[1][1]);
                M[0][0] /= det;
                M[1][1] /= det;
                M[0][1] /= -det;
                M[1][0] /= -det;
                return true;
            }

            void PseudoInverseMatrix2x2(Matrix2x2 &in) {
                Matrix2x2 M,inT;
                TransposeMatrix2x2(in,inT);
                MultiplyMatrix2x2(inT,in,M);
                InverseMatrix2x2(M);
                MultiplyMatrix2x2(M,inT,in);
            }
        };   // end of the declaration of the class "Hamiltonian". =====================

        // Hamiltonian
        Hamiltonian *hamiltonian;

    public:
        // Constructor
        AnisotropicTensorDistance(T *_data, int _width, int _height, int _depth, T* mask, T *tensor, T _dx = T(1), T _dy = T(1), T _dz = T(1)) :
            PradosSchemesForFastMarching_3D<T>(_data,_width,_height,_depth), dx(_dx), dy(_dy), dz(_dz) {

            // Allocation and initialization of the Hamiltonian from the tensor field
            hamiltonian = new Hamiltonian[this->size];
            for (int x=0;x<this->width;x++) for (int y=0;y<this->height;y++) for(int z=0;z<this->depth;z++) {
                const int n = this->_offset(x,y,z);
                if (mask[n]) {
                    Matrix3x3 B = { { tensor[n], tensor[n+this->size], tensor[n+2*this->size] },
                                    { tensor[n+this->size], tensor[n+3*this->size], tensor[n+4*this->size] },
                                    { tensor[n+2*this->size], tensor[n+4*this->size], tensor[n+5*this->size] }};
                                    if (!hamiltonian[n].InitFromB(B)) {
                                        // std::cerr << "Null tensor was found at ("<< x << ", " << y << ", " << z << "), discarding from mask..." << std::endl;
                                        mask[n] = 0;
                                        this->AddForbiddenPoint(x,y,z);
                                    }
                } else this->AddForbiddenPoint(x,y,z);
            }
        }


        // Destructor
        virtual ~AnisotropicTensorDistance() {
            delete [] hamiltonian;
        }

        virtual void setOptDynamics(int x, int y, int z, T component1, T component2, T component3) const {
            const int n = this->_offset(x,y,z);
            const Matrix3x3 &iB = hamiltonian[n].invB;
            const T &norm = std::sqrt(iB[0][0]*component1*component1 + iB[1][1]*component2*component2 + iB[2][2]*component3*component3 + 2*iB[0][1]*component1*component2 + 2*iB[0][2]*component1*component3 + 2*iB[1][2]*component2*component3);
            this->OptDynamics[n][0] = component1/norm;
            this->OptDynamics[n][1] = component2/norm;
            this->OptDynamics[n][2] = component3/norm;
        }

    protected:
        virtual bool eqSolverOnPart_with_s1s2s3_nonNull(
            const T U1, const T U2, const T U3,             // Values of the solution at the considered neigborhood voxels,
            const int s1, const int s2, const int s3,       // signs associated to the considered sector,
            const int x, const int y, const int z,          // coordinates of the considered voxel,
            T &Root,                                        // solution.
            T &optDymamics1,    T &optDymamics2,    T &optDymamics3 // optimal dynamic associated to the solution.
            ) const {

            if (U1==this->big || U2==this->big || U3==this->big) return false; // of course there is no such solution :-)

            // Solving of the equation
            // p_t^T B p_t - 1 = 0

            T t1=this->big, t2=this->big;

            const Matrix3x3 &B = hamiltonian[this->_offset(x,y,z)].B;

            if (basicAnisoEikonalEq_3D(
                B[0][0], B[0][1], B[0][2],
                         B[1][1], B[1][2],
                                  B[2][2],
                U1,U2,U3,
                s1*dx,s2*dy,s3*dz,
                t1, t2)) {

                // test of p_t1:            ($[p_t]_i = [ t-u(x+sihiei)] / (-sihi)$.)
                T  // description of p_t1
                    p_t1_1 =  ( t1-U1 ) / (-s1*dx),
                    p_t1_2 =  ( t1-U2 ) / (-s2*dy),
                    p_t1_3 =  ( t1-U3 ) / (-s3*dz);

                // with this equation, the optimal dymamics $f(x,a_p)$ associated to a vector $p$  is $B*p$.
                // So the optimal dynamics associated to $t1$ is $B*p_t1$:
                optDymamics1 =
                    // Bp_t1_1 =
                                B[0][0]*p_t1_1 + B[0][1]*p_t1_2 + B[0][2]*p_t1_3;
                optDymamics2 =
                    // Bp_t1_2 =
                                B[0][1]*p_t1_1 + B[1][1]*p_t1_2 + B[1][2]*p_t1_3;
                optDymamics3 =
                    // Bp_t1_3 =
                                B[0][2]*p_t1_1 + B[1][2]*p_t1_2 + B[2][2]*p_t1_3;

                // if for all i, [B*p_t1]_i < 0 then t1 is the solution
                if ( (optDymamics1*s1 < 0) && (optDymamics2*s2 < 0) && (optDymamics3*s3 < 0) ) {
                    Root = t1;
                    return true;
                }

                // test of p_t2:
                T
                    p_t2_1 =  ( t2-U1 ) / (-s1*dx),
                    p_t2_2 =  ( t2-U2 ) / (-s2*dy),
                    p_t2_3 =  ( t2-U3 ) / (-s3*dz);

                // The optimal dynamics associated to $t2$ is $B*p_t2$:
                optDymamics1 =
                    //  Bp_t2_1 =
                                B[0][0]*p_t2_1 + B[0][1]*p_t2_2 + B[0][2]*p_t2_3;
                optDymamics2 =
                    //  Bp_t2_2 =
                                B[0][1]*p_t2_1 + B[1][1]*p_t2_2 + B[1][2]*p_t2_3;
                optDymamics3 =
                    //  Bp_t2_3 =
                                B[0][2]*p_t2_1 + B[1][2]*p_t2_2 + B[2][2]*p_t2_3;

                // if for all i, [B*p_t2]_i < 0 then t2 is the solution
                if ( (optDymamics1*s1 < 0) && (optDymamics2*s2 < 0) && (optDymamics3*s3 < 0) ) {
                    Root = t2;
                    return true;
                }
            }

            // There is no solution
            return  false;
        }


        virtual bool eqSolverOnPart_withOne_si_Null(
            const T U1, const T U2, const T U3,
            const int s1, const int s2, const int s3,
            const int x, const int y, const int z,
            const int indice_si_EqualZero,
            T &Root,
            T &optDymamics1,    T &optDymamics2,    T &optDymamics3
            ) const {

            const int n = this->_offset(x,y,z);

            if (indice_si_EqualZero==1) {
                if (U2==this->big || U3==this->big) return false; // of course there is no such solution :-)
                const Matrix2x2 &invinvB = hamiltonian[n].invinvB1;
                return solveInDim2(invinvB[0][0], invinvB[0][1], invinvB[1][1], U2, U3, s2, s3, dy, dz, Root,optDymamics2,optDymamics3);
            }

            if (indice_si_EqualZero==2) {
                if (U1==this->big || U3==this->big) return false; // of course there is no such solution :-)
                const Matrix2x2 &invinvB = hamiltonian[n].invinvB2;
                return solveInDim2(invinvB[0][0], invinvB[0][1], invinvB[1][1], U1, U3, s1, s3, dx, dz, Root,optDymamics1,optDymamics3);
            }

            if (indice_si_EqualZero==3) {
                if (U1==this->big || U2==this->big) return false; // of course there is no such solution :-)
                const Matrix2x2 &invinvB = hamiltonian[n].invinvB3;
                return solveInDim2( invinvB[0][0], invinvB[0][1], invinvB[1][1], U1, U2, s1, s2, dx, dy, Root,optDymamics1,optDymamics2);
            }

            // In any case, indice_si_EqualZero must be equal to 1, 2 or 3 !!!!!
            assert(false);
            return  false;
        }


        virtual bool eqSolverOnPart_withTwo_si_Null(
            const T U1, const T U2, const T U3,
            const int s1, const int s2, const int s3,
            const int x, const int y, const int z,
            const int indice_si_DiffZero,
            T &Root,
            T &optDymamics1,    T &optDymamics2,    T &optDymamics3
            ) const {

            const int n = this->_offset(x,y,z);

            if (indice_si_DiffZero==1) {
                if (U1==this->big) return false; // of course there is no such solution :-)
                return solveInDim1(hamiltonian[n].invsqinvA1,U1,s1,dx,Root,optDymamics1);
            }

            if (indice_si_DiffZero==2) {
                if (U2==this->big) return false; // of course there is no such solution :-)
                return solveInDim1(hamiltonian[n].invsqinvA2,U2,s2,dy,Root,optDymamics2);
            }

            if (indice_si_DiffZero==3) {
                if (U3==this->big) return false; // of course there is no such solution :-)
                return solveInDim1(hamiltonian[n].invsqinvA3,U3,s3,dz,Root,optDymamics3);
            }

            // In any case, indice_si_DiffZero must be equal to 1, 2 or 3 !!!!!
            assert(false);
            return  false;
        }


        bool solveInDim2(
            const T C00, const T C01,
                         const T C11,
            const T U1, const T U2,
            const int s1, const int s2,
            const T dx1, const T dx2,
            T &root,
            T &optDymamics1,    T &optDymamics2) const {

            // Solving of the equation
            // p_t^T C p_t - 1 = 0

            T t1=this->big, t2=this->big;

            if (basicAnisoEikonalEq_2D(C00, C01, C11, U1, U2, s1*dx1, s2*dx2, t1, t2)) {

                // test of p_t1:            ($[p_t]_i = [ t-u(x+sihiei)] / (-sihi)$.)
                T  // description of p_t1
                    p_t1_1 =  ( t1-U1 ) / (-s1*dx1),
                    p_t1_2 =  ( t1-U2 ) / (-s2*dx2);

                // The optimal dynamics associated to $t1$ is $C*p_t1$:
                optDymamics1 = // Cp_t1_1 =
                                            C00*p_t1_1 + C01*p_t1_2;
                optDymamics2 = // Cp_t1_2 =
                                            C01*p_t1_1 + C11*p_t1_2;

                // if for all i, [C*p_t1]_i < 0 then t1 is the solution
                if ( (optDymamics1*s1 < 0) && (optDymamics2*s2 < 0) ) {
                    root = t1;
                    return true;
                }

                // test of p_t2:
                T   p_t2_1 =  ( t2-U1 ) / (-s1*dx1),
                    p_t2_2 =  ( t2-U2 ) / (-s2*dx2);

                // The optimal dynamics associated to $t2$ is $C*p_t2$:
                optDymamics1 = // Cp_t2_1 =
                                            C00*p_t2_1 + C01*p_t2_2;
                optDymamics2 = // Cp_t2_2 =
                                            C01*p_t2_1 + C11*p_t2_2;

                // if for all i, [B*p_t2]_i < 0 then t2 is the solution
                if ( (optDymamics1*s1 < 0) && (optDymamics2*s2 < 0) ) {
                    root = t2;
                    return true;
                }
            }

            // There is no solution
            return  false;
        }


        bool solveInDim1(
            const T C00,
            const T U1,
            const int s1,
            const T dx1,
            T &root,
            T &optDymamics1) const {

            // The 1D case is particularly simple
            if (C00>0) {
                T sqrtC00 = (T)std::sqrt(C00);
                root = U1 + dx1 / sqrtC00;
                optDymamics1 = -s1 * sqrtC00;
                return true;
            }

            return false;
        }


        //////////////////////////////////////////////////////////////////////////
        // Inversion of the eikonal equation in a given simplex

        bool basicAnisoEikonalEq_2D(
            T c11, T c12,               // Matrix C
                T c22,
            T u1, T u2,                 // Some values of U
            T dx1, T dx2,               // Signed mesh size
            T &sol_max,  T &sol_min) const {

            const T dx1_2 = dx1*dx1;
            const T dx2_2 = dx2*dx2;

            const T d1 = dx2_2*c11;
            const T d2 = dx1_2*c22;
            const T d12 = 2*dx1*dx2*c12;

            const T a = d1 + d2 + d12;
            const T b = - T(2)*d1*u1 - T(2)*d2*u2 - d12*(u1+u2);
            const T c = d1*u1*u1 + d2*u2*u2 + d12*u1*u2 - dx1_2*dx2_2;

            return this->_SolveTrinome(a,b,c,sol_max,sol_min);
        }


        bool basicAnisoEikonalEq_3D(
            T c11, T c12, T c13,        // Matrix C
                T c22, T c23,
                        T c33,
            T u1, T u2, T u3,           // Some values of U
            T dx1, T dx2, T dx3,        // Signed mesh size
            T &sol_max,  T &sol_min) const {

            const T dx1_2 = dx1*dx1;
            const T dx2_2 = dx2*dx2;
            const T dx3_2 = dx3*dx3;

            const T d1 = dx2_2*dx3_2*c11;
            const T d2 = dx1_2*dx3_2*c22;
            const T d3 = dx1_2*dx2_2*c33;

            const T d12 = 2*dx1*dx2*dx3_2*c12;
            const T d13 = 2*dx1*dx2_2*dx3*c13;
            const T d23 = 2*dx1_2*dx2*dx3*c23;

            const T a = d1 + d2 + d3 + d12 + d13 + d23;
            const T b = - 2*d1*u1 - 2*d2*u2 - 2*d3*u3 - d12*(u1+u2) - d13*(u1+u3) - d23*(u2+u3);
            const T c = d1*u1*u1 + d2*u2*u2 + d3*u3*u3 + d12*u1*u2 + d13*u1*u3 + d23*u2*u3 - dx1_2*dx2_2*dx3_2;

            return this->_SolveTrinome(a,b,c,sol_max,sol_min);
        }

    }; // End of the class AnisotropicTensorDistance.

} // End of the namespace "FastLevelSet".

#endif

