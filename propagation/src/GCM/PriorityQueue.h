/*------------------------------------------------------------------------------------------------------
  
  File        : PriorityQueue.h   (GCM Library)

  Authors     : Emmanuel Prados (UCLA/INRIA), Christophe Lenglet (INRIA), Jean-Philippe Pons (INRIA)
  
  Description : This file contains severals classes allowing to manage with the  priority queue. 

  --------------

  License     : This software is governed by the CeCILL-C license under French law and abiding by the rules of distribution of free software. 

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
    
  
  ------------------------------------------------------------------------------------------------


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



#ifndef FASTLEVELSETPRIORITYQUEUE_H
#define FASTLEVELSETPRIORITYQUEUE_H

#include <vector>

namespace FastLevelSet {

    template <class value_type>
    class PriorityQueue;
    

    // Class for a node of the binary tree implementing the queue
    template <class value_type>
        class PriorityQueueNode {
        value_type value;
        PriorityQueueNode *father, *first_son, *next_son;
        PriorityQueueNode(const value_type &x) : value(x), father(0), first_son(0), next_son(0) {}
        friend class PriorityQueue<value_type>;
    };


    // Class for a priority queue with the ability to increase the priority of items in the queue
    template <class value_type>
    class PriorityQueue {

      public:
        // Constructs an empty priority queue
        PriorityQueue() : subtree_tab(5) { root = 0; }
        
        
        // Destructs the priority queue
        ~PriorityQueue() { clear(); }
        
        
        // Returns true if the priority queue is empty, false otherwise
        bool empty() const { return (root==0); }
    
    
        // Returns a constant reference to the element in the queue with the highest priority
        const value_type & top() const { return root->value; }
        
        
        // Adds x to the queue
        PriorityQueueNode<value_type> *push(const value_type &x) {
            
            PriorityQueueNode<value_type> *new_node = new PriorityQueueNode<value_type>(x);
            
            if (empty())
                root = new_node;
            else
                _insert_node(root, new_node);

            return new_node;
        }
        
        
        // Removes the item with the highest priority from the queue
        void pop() {

            if (empty()) return;
            
            PriorityQueueNode<value_type> *old_root = root;
            if (root->first_son == 0)
                root = 0;
            else
                root = _combine_sons(root->first_son);
            
            delete old_root;
        }

        
        // Deletes all elements from the priority queue
        void clear() {
            _clear(root);
            root = 0;
        }
        
        
        // Increases the priority of an item of the queue
        void increase(PriorityQueueNode<value_type> *node, const value_type &x) {
            if (node->value < x) return; // The priority can only increase
            
            node->value = x;
            if (node != root) {
                if (node->next_son != 0)
                    node->next_son->father = node->father;
                if (node->father->first_son == node)
                    node->father->first_son = node->next_son;
                else
                    node->father->next_son = node->next_son;
                
                node->next_son = 0;
                _insert_node(root, node);
            }
        }
        
        
      private:
        // Root of the tree
        PriorityQueueNode<value_type> *root;

        // Subtree array
        std::vector<PriorityQueueNode<value_type> *> subtree_tab;
        
        // Recursively deletes the nodes of the tree
        void _clear(PriorityQueueNode<value_type> *node) const {
            if (node) {
                _clear(node->first_son);
                _clear(node->next_son);
                delete node;
            }
        }
        
    
        void _insert_node(PriorityQueueNode<value_type> *&first, PriorityQueueNode<value_type> *second) const {
            
            if (second == 0) return;
            
            if (second->value < first->value) {
                second->father = first->father;
                first->father = second;
                first->next_son = second->first_son;
                if (first->next_son != 0)
                    first->next_son->father = first;
                second->first_son = first;
                first = second;
            }
            else {
                second->father = first;
                first->next_son = second->next_son;
                if (first->next_son != 0)
                    first->next_son->father = first;
                second->next_son = first->first_son;
                if (second->next_son != 0)
                    second->next_son->father = second;
                first->first_son = second;
            }
        }
        
        
        PriorityQueueNode<value_type> *_combine_sons(PriorityQueueNode<value_type> *first) {
            
            if (first->next_son == 0) return first;
            
            int n = 0;
            for(;first!=0;n++) {
                if (n == int(subtree_tab.size()))
                    subtree_tab.resize(n*2);
                subtree_tab[n] = first;
                first->father->next_son = 0;
                first = first->next_son;
            }
            if (n == int(subtree_tab.size()))
                subtree_tab.resize(n+1);
            subtree_tab[n] = 0;
            
            int i = 0;
            for(;i+1<n;i+=2)
                _insert_node(subtree_tab[i], subtree_tab[i+1]);
            
            int j = i-2; 
            if (j==n-3)
                _insert_node(subtree_tab[j], subtree_tab[j+2]);
          
            for(;j>=2;j-=2)
                _insert_node(subtree_tab[j-2], subtree_tab[j]);
            return subtree_tab[0];
        }
    };
}

#endif
