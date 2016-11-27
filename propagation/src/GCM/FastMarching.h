/*------------------------------------------------------------------------------------------------------
  
  File        : FastMarching.h   (GCM Library)

  Description : This class contains all the basic tools for working the fast marching algorithm.

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


#ifndef FASTLEVELSETFASTMARCHING_H
#define FASTLEVELSETFASTMARCHING_H

#include <algorithm>
#include <limits>
#include "Globals.h"
#include "PriorityQueue.h"

namespace FastLevelSet {
    typedef enum { eAlive=0, eTrial=1, eFar=2, eForbidden=3 } eState;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // English: Class containing the coordinates of a point and the value of the distance function at this point. 
    // French: Classe contenant les coordonnees d'un point et la valeur de la fonction distance en ce point
    template <typename T>
    struct stCoordVal : public stCoord {
        T val;
        stCoordVal(int _x, int _y, int _z, int _n, T _val) : stCoord(_x,_y,_z,_n), val(_val) {}
        stCoordVal(const stCoord &coord, T _val) : stCoord(coord), val(_val) {}

        bool operator < (const stCoordVal &a) const {
            return (val < a.val);
        }
    };


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // En: Class for the fast marching    
    // Fr: Classe pour le fast marching    
    template <typename T = float, int sign = +1>
    class FastMarching {
      private:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Types
        // Fr: Types
        typedef PriorityQueue< stCoordVal<T> > CoordQueue;
        typedef PriorityQueueNode< stCoordVal<T> > QueueNode;
          
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Members
        // Fr: Membres
        unsigned char *state;   // En: state of the points // Fr: Etat des points
        QueueNode **tabnode;    // En: Tab of nodes of the queue : needed for modifying the priority of a ponit of the queue // Fr. :Tableau des noeuds de la queue : necessaire pour modifier la priorite d'un point de la queue
        CoordQueue trial_queue;  // En: Queue of the trial points // fr: Queue des points trial
		CoordList alive_list;   // En: list of the alive points // Fr: Liste des points alive
        T limit;                // En: stopping test // Fr: Condition d'arret du fast marching  
      protected:
        T *data;                // En: Tab containing the data // Fr: Tableau des donnees
        int width,height,depth,size; // En: Dimension of the data // Fr: Dimension des donnees
        T big;                 // En: Maximal value authorized for the type T // Fr: Valeur maximale permise par le type T

      public:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Constructor / destructor
        // Fr: Constructeur / destructeur
		FastMarching(T *_data, int _width, int _height, int _depth = 1) : data(_data), width(_width), height(_height), depth(_depth), size(width*height*depth), big(std::numeric_limits<T>::max()) {
            state = new unsigned char[size];
            tabnode = new QueueNode *[size];
            Init();
        }
        virtual ~FastMarching() {
            delete state;
            delete tabnode;
        }

		////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Iterators on the list of alive points
        // Fr: Iterateurs sur la liste de points alive
		const_CoordListIterator begin() const { return alive_list.begin(); }
        const_CoordListIterator end() const { return alive_list.end(); }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Access to the state of a point 
        // Fr: Acces a l'etat d'un point
        eState GetState(const int x, const int y, const int z = 0) const { return GetState(_offset(x,y,z)); }
        eState GetState(const int n) const { return eState(state[n]); }


	////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Initialization
        // Fr: Initialisation
		void Init() {
			trial_queue.clear();
			alive_list.clear();
            memset(state,eFar,size*sizeof(unsigned char));
            memset(tabnode,0,size*sizeof(QueueNode*));
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Addition of a trial point in the queue
        // Fr: Ajout d'un point trial a la queue
	void AddTrialPoint(const int x, const int y, const int z = 0) { AddTrialPoint(x,y,z,_offset(x,y,z)); }
	void AddTrialPoint(const int x, const int y, const int z, const int n) {
		    if (state[n] == eFar) {
		       	state[n] = eTrial;
		       	tabnode[n] = trial_queue.push(stCoordVal<T>(x,y,z,n,sign*data[n]));
		    }
        }


	////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En : Addition of an alive point to the list 
        // Fr : Ajout d'un point alive a la liste
	void AddAlivePoint(const int x, const int y, const int z = 0) { AddAlivePoint(x,y,z,_offset(x,y,z)); }
	void AddAlivePoint(const int x, const int y, const int z, const int n) {
	       	if (state[n] == eFar) {
		       	state[n] = eAlive;
		       	alive_list.push_back(stCoord(x,y,z,n));
	       	}
        }

        
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Addition of a forbidden point to the list 
        // Fr: Ajout d'un point interdit
	void AddForbiddenPoint(const int x, const int y, const int z = 0) { AddForbiddenPoint(_offset(x,y,z)); }
        void AddForbiddenPoint(int n) {
            if (state[n] == eFar) {
                state[n] = eForbidden;
            }
        }


	////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Generation of the trial points from the alive points 
        // Fr: Generation des points trial a partir des points alive
       	void InitTrialFromAlive() {
	       	// En: We read the list of alive points and we add their neibourghood points to the queue of the trial points
	       	// Fr: On lit la liste de points alive et on ajoute les voisins a la queue des points trial
	       	for (const_CoordListIterator vox=begin();vox!=end();vox++) {
	       		const int x = vox->x;
			const int y = vox->y;
			const int z = vox->z;
			if (x>0) AddTrialPoint(x-1,y,z);
			if (x<width-1) AddTrialPoint(x+1,y,z);
			if (y>0) AddTrialPoint(x,y-1,z);
			if (y<height-1) AddTrialPoint(x,y+1,z);
			if (z>0) AddTrialPoint(x,y,z-1);
			if (z<depth-1) AddTrialPoint(x,y,z+1);
	       	}
       	}


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Launching of the fast marching
        // Fr: Lancement du fast marching
       	void Run() {
	       	Run( (sign>0) ? big : -big );
       	}

       	void Run(T _limit) {
            // En: We memorize the spread of the fast marching
            // Fr: On memorise l'etendue du fast marching
            limit = _limit;

            // En: As much as there are trial points in the queue...
            // Fr: Tant qu'il y a des points trial dans la queue
            while(!trial_queue.empty()) {
                    
                // En: We are looking for the trial point which has the smallest value
                // Fr: On cherche le point trial de valeur la plus faible
                stCoordVal<T> pt = trial_queue.top();
                                
                // En: If we reach the limit
                // Fr: Si on a atteint la limite
                if (pt.val>=sign*limit) {
                    // En: we empty the queue and we put all points as far 
                    // Fr: On vide la queue et on marque tous ses points comme far
                    while (!trial_queue.empty()) {
                        pt = trial_queue.top();
                        trial_queue.pop();
                        state[pt.n] = eFar;
                        tabnode[pt.n] = 0;
                        data[pt.n] = limit;
                    } 
                    return;
                }
                
                // En: We remove the point of the queue and we add it to the list of alive points
                // Fr: On enleve le point de la queue et on l'ajoute a la liste des points alive
                trial_queue.pop();
                const int n = pt.n;
				state[n] = eAlive;
                tabnode[n] = 0;
				alive_list.push_back(pt);
                
                // En: we determine the neibourghood points and we update them 
                // fr: On determine les voisins du point et on les met a jour
				const int x = pt.x;
                const int y = pt.y;
                const int z = pt.z;
                if (x>0) __UpdatePoint(x-1,y,z);
                if (x<width-1) __UpdatePoint(x+1,y,z);
                if (y>0) __UpdatePoint(x,y-1,z);
                if (y<height-1) __UpdatePoint(x,y+1,z);
                if (z>0) __UpdatePoint(x,y,z-1);
                if (z<depth-1) __UpdatePoint(x,y,z+1);
            }
        }


    private:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Update of a point
        // Fr: Mise a jour d'un point
       void __UpdatePoint(const int x, const int y, const int z) {
			const int n = _offset(x,y,z);
            const eState st = (eState)state[n];
            if (st == eFar) {
                // En: If it was a far point, we put it as trial,  we compute its value and we add it to the queue
                // Fr: Si c'etait un point far, on le rend trial, on calcule sa valeur et on l'ajoute a la queue
                state[n] = eTrial;
                const T val = _UpdateValue(x,y,z);
                data[n] = sign*val;
				tabnode[n] = trial_queue.push(stCoordVal<T>(x,y,z,n,val));
            }
            else if (st == eTrial) {
                // En: If it was a far point, we recompute its value
                // En: if the new value is inferior, we accept it and we increase the point in the queue
                // Fr: Si c'etait deja un point trial, on recalcule sa valeur
                // Fr: Si la nouvelle valeur est inferieure, on l'accepte et on fait remonter le point dans la queue
                const T val = _UpdateValue(x,y,z);
                if (val<sign*data[n]) {
                    data[n] = sign*val;
                    trial_queue.increase(tabnode[n],stCoordVal<T>(x,y,z,n,val));
                }
            }
        }

	protected:
      	////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Acces to a  point
        // Fr: Acces a un point
        int _offset(const int x, const int y, const int z) const { return x+width*(y+height*z); }
        T _GetValue(const int x, const int y, const int z) const { return _GetValue(_offset(x,y,z)); }
        T _GetValue(const int n) const {
            if (state[n] <= eTrial) return sign*data[n];
            return big;
		}

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Fr: Computation of the value of a point : have to be defined in the inherited class
        // En: Calcul de la valeur d'un point : a definir dans les classes derivees
        virtual T _UpdateValue(const int x, const int y, const int z) const = 0;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: resolution of a second degree trinomial equation
        // Fr: Resolution d'un trinome du second degre
        bool _SolveTrinome(const T a, const T b, const T c, T &sol_max, T &sol_min) const {
            const T delta = b*b - 4*a*c;
            if (delta < 0) return false;
            const T sqrtDelta = std::sqrt(delta);
            sol_max = (- b + sqrtDelta) / a / 2;
            sol_min = (- b - sqrtDelta) / a / 2;
            return true;
        }
    };

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // En: 2D Iconale Equation 
    // Fr: Equation iconale en 2D
    template <typename T = float, int sign = +1>
    class Eikonal2D : public FastMarching<T,sign> {
            
      public:
        Eikonal2D(T *_data, int _width, int _height) :  FastMarching<T,sign>(_data,_width,_height) {}
        virtual ~Eikonal2D() {}
        
      protected:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////     
        // En: Update of a point
        // Fr: Mise a jour d'un point
        virtual T _UpdateValue(const int x, const int y, const int z) const {
   
            // En: we take the minimum of each couple of neighbourhood
            // Fr: On prend le minimum de chaque paire de voisins
            const T A = (x==0) ? this->_GetValue(x+1,y,z) : (x==this->width-1) ? this->_GetValue(x-1,y,z) : std::min( this->_GetValue(x+1,y,z), this->_GetValue(x-1,y,z) );
            const T B = (y==0) ? this->_GetValue(x,y+1,z) : (y==this->height-1) ? this->_GetValue(x,y-1,z) : std::min( this->_GetValue(x,y+1,z), this->_GetValue(x,y-1,z) );

            return _Solve(A,B);
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Resolution of the trinomial equation
        // Fr: Resolution du trinome
        T _Solve(T A, T B) const {

            // En: we get rid of the trival cases 
            // Fr: On se debarasse des cas triviaux
            if (A == this->big) return B+1;
            if (B == this->big) return A+1;

            // En: we reorder the values in order to have  B>=A
            // Fr: On reordonne les valeurs pour avoir B>=A
            if (A>B) std::swap(A,B);
           
            //En: We assume that  u>=B : we solve the trinomial equation. If we have u>=B, it is ok
            //Fr: On supppose u>=B : trinome. Si on a bien u>=B, on a gagne
            T sol_max, sol_min;
            if (B<this->big && _SolveTrinome(2, -2*(A+B), A*A+B*B-1, sol_max, sol_min) && sol_max+EPS>=B) return sol_max;

            // En: We assume A<=u<B
            // Fr: On suppose A<=u<B
            return A+1;
        }
    };

    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // En: 3D Iconale Equation 
    // Fr: Equation iconale en 3D
    template <typename T = float, int sign = +1>
    class Eikonal3D : public FastMarching<T,sign> {
    
    public:
        Eikonal3D(T *_data, int _width, int _height, int _depth) :  FastMarching<T,sign>(_data,_width,_height,_depth) {}
        virtual ~Eikonal3D() {}

    protected:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Update of a point
        // Fr: Mise a jour d'un point
        virtual T _UpdateValue(const int x, const int y, const int z) const {

            // En: we take the minimum of each couple of neighbourhood
            // Fr: On prend le minimum de chaque paire de voisins
            const T A = (x==0) ? this->_GetValue(x+1,y,z) : (x==this->width-1) ? this->_GetValue(x-1,y,z) : std::min( this->_GetValue(x+1,y,z), this->_GetValue(x-1,y,z) );
            const T B = (y==0) ? this->_GetValue(x,y+1,z) : (y==this->height-1) ? this->_GetValue(x,y-1,z) : std::min( this->_GetValue(x,y+1,z), this->_GetValue(x,y-1,z) );
            const T C = (z==0) ? this->_GetValue(x,y,z+1) : (z==this->depth-1) ? this->_GetValue(x,y,z-1) : std::min( this->_GetValue(x,y,z+1), this->_GetValue(x,y,z-1) );

            return _Solve(A,B,C);
        }


        ////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // En: Resolution of the trinomial equation
        // Fr: Resolution du trinome
        T _Solve(T A, T B, T C) const {

            // En: we reorder the values in order to have  C>=B>=A
            // Fr: On reordonne les valeurs pour avoir C>=B>=A
            if (A>B) std::swap(A,B);
            if (B>C) std::swap(B,C);
            if (A>B) std::swap(A,B);

            // En: We assume sol>=C : first trinomial equation. If we have sol>=C, it is ok
            // Fr: On suppose sol>=C : premier trinome. Si on a bien sol>=C, on a gagne
            T sol_max, sol_min;
            if (C<this->big && _SolveTrinome(3, -2*(A+B+C), A*A+B*B+C*C-1, sol_max, sol_min) && sol_max+EPS>=C) return sol_max;

            // En: We assume B<=sol<C : second trinomial equation. If we have sol>=B, it is ok
            // Fr: On supppose B<=sol<C : deuxieme trinome. Si on a bien sol>=B, on a gagne
            if (B<this->big && _SolveTrinome(2, -2*(A+B), A*A+B*B-1, sol_max, sol_min) && sol_max+EPS>=B) return sol_max;
                
            // En: We assume A<=sol<B
            // Fr: On suppose A<=sol<B
            return A+1;
        }
    };
}

#endif
