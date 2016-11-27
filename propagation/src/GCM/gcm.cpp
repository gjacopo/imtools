/*------------------------------------------------------------------------------------------------------
  
  File        : gcm.cpp

  Description : Example of use of the GCM library

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

  Associated publications  : This C++ code corresponds to the implementation of the algorithm presented in the following articles:  
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; Control Theory and Fast Marching Methods for Brain Connectivity Mapping; INRIA Research Report 5845 -- UCLA Computer Science Department Technical Report 060004, February 2006.
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; Control Theory and Fast Marching Methods for Brain Connectivity Mapping; Proc. IEEE Computer Society Conference on Computer Vision and Pattern Recognition, New York, NY, I: 1076-1083, June 17-22, 2006.  
  - For more references, we refer to the official web page of the GCM Library and to authors' web pages.

  Please, if you use the GCM library in your work, make sure you will include the reference to the work of the authors in your publications.
  
  --------------
  
  Technical Comments and Detailled description:
   - Goal: example of use of the GCM library.
   - "CImg" library: to load and manage (Diffusion Tensor) Images, this example uses the "CImg" library: http://cimg.sourceforge.net/ 
     It has been compiled and tested with the version 1.0.7 of the "CImg" library.
   Note: even if this example file requires the "CImg" library, the GCM library is completely independent on any image library, in particular it is independent on the CImg Library.
   - Tasks:
   * loading of the data (image, starting point, optional mask)
   * initialization of the objects
   * Computation the distance function/optimal dynamics/confidence statistics within the mask
   * Saving distance function/optimal dynamics/confidence statistics
  
  ----------------------------------------------------------------------------------------------------*/


#include <iostream>
#include <cfloat>
#include <climits>
#include "CImg.h"
#include "AnisotropicTensorDistance.h"
#include "AnisotropicTensorDistanceConfidence.h"
using namespace cimg_library;
using namespace std;
using namespace FastLevelSet;

typedef float T;

int main(int argc, char **argv) {

    // Parameters
    cimg_usage("Geodesic Connectivity Mapping\n");

    const char* tensors_name = cimg_option("-dti",(const char *)0,"DTI");
    const char* origin = cimg_option("-orig","0,0,0","Origin");
    const char* mask_name = cimg_option("-mask",(const char *)0,"Mask");
    const T alpha = cimg_option("-alpha",0.0,"Metric power");
    const char* distance_name = cimg_option("-dist","","Distance output");
    const char* odyn_name = cimg_option("-odyn","","Optimal dynamics output");
    const char* cmean_name = cimg_option("-mean","","Confidence mean output");
    const char* cstd_name = cimg_option("-std","","Confidence std output");
    const char* cmin_name = cimg_option("-min","","Confidence min output");

    if (cimg_option("-h",(const char *)0,0)) return 0;
    if (argc<4) {
        cerr << "Not enough arguments" << endl;
        return 1;
    }

    // Load the tensor field
    cout << "* Loading tensor field from '" << tensors_name << "'... "; cout.flush();
    CImg<T> tensor(tensors_name);
    cout << tensor.width << " x " << tensor.height << " x " << tensor.depth << " x " << tensor.dim << " "; cout.flush();
    cout << "OK" << endl;

    // Load the mask
    CImg<T> mask;
    if (mask_name) {
        cout << "* Loading mask from '" << mask_name << "'... "; cout.flush();
        mask = CImg<T>::load(mask_name);
        cout << mask.width << " x " << mask.height << " x " << mask.depth << " "; cout.flush();
    }
    cout << "OK" << endl;

    // Tensors power
    CImg<T> tensor_power = tensor;

    // tensor_power already contains tensor and if alpha==0 the code in
    // AnisotropicTensorDistanceConfidence.h (l. 62) will explicitely compute the Euclidean norm of the dynamics
    if (alpha != 0.0 && alpha != 1.0) {
        cout << "* Computing tensors power for confidence estimation ... "; cout.flush();

	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	// Note: For executing this part, you need to link with LAPACK in compilation process.
	// See http://www.netlib.org/lapack/
	// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        CImg<T> val(3,1,1,1), valM(3,3,1,1), vec(3,3,1,1);
        valM.fill(0.0);
		vec.fill(0.0);

        cimg_mapXYZ(tensor_power,x,y,z) {
            const CImg<T> D = tensor_power.get_tensor(x,y,z);
			if (D(0,0)>std::numeric_limits<T>::epsilon()) {
				D.symeigen(val,vec);
				valM(0,0) = std::pow(val(0),alpha);
				valM(1,1) = std::pow(val(1),alpha);
				valM(2,2) = std::pow(val(2),alpha);
				tensor_power.set_tensor(vec.get_transpose()*valM*vec,x,y,z);
			}
		}
        cout << "OK" << endl;
    }


    // Fast marching initialization
    cout << "* Initializing... "; cout.flush();
    CImg<T> dist(tensor.width,tensor.height,tensor.depth);
    dist.fill(FLT_MAX);
    AnisotropicTensorDistanceConfidence<T> march(dist.data,dist.width,dist.height,dist.depth,mask.data,tensor.data,tensor_power.data,alpha);
    cout << "OK" << endl;


    int x0,y0,z0;
    sscanf(origin,"%d,%d,%d",&x0,&y0,&z0);
    dist(x0,y0,z0) = 0;
    march.AddTrialPoint(x0,y0,z0);


    // Compute the distance function, the optimal dynamics and confidence statistics, within the mask
    cout << "* Computing distance... "; cout.flush();
    march.Run();
    cout << "OK" << endl;

    // Save the result
    if (distance_name != "") {
        cout << "* Saving distance to '" << distance_name << "'... "; cout.flush();
        dist.save(distance_name);
        cout << "OK" << endl;
    }

    // Extract optimal dynamics   (the optimal dynamics have been already computed; the result of the optimal dynamics is stacked in some variables of the "march" object)
    CImg<T> optDynamics(tensor.width,tensor.height,tensor.depth,3);
    cimg_mapXYZ(optDynamics,x,y,z) {
        optDynamics(x,y,z,0) = march.getOptDynamics(x,y,z,0);
        optDynamics(x,y,z,1) = march.getOptDynamics(x,y,z,1);
        optDynamics(x,y,z,2) = march.getOptDynamics(x,y,z,2);
    }

    if (odyn_name != "") {
        cout << "* Saving optimal dynamics to '" << odyn_name << "'... "; cout.flush();
        optDynamics.save(odyn_name);
        cout << "OK" << endl;
    }

    // Get confidence   (as the optimal dynamics, the confidence has been already computed!)
    CImg<T> confidenceMean(tensor.width,tensor.height,tensor.depth);
    CImg<T> confidenceStd(tensor.width,tensor.height,tensor.depth);
    CImg<T> confidenceMin(tensor.width,tensor.height,tensor.depth);

    cimg_mapXYZ(confidenceMean,x,y,z) {
        if (x==x0 && y==y0 && z==z0) {
            confidenceMean(x,y,z) = 0.0;
            confidenceStd(x,y,z) = 0.0;
            confidenceMin(x,y,z) = 0.0;
        } else {
            confidenceMean(x,y,z) = march.getConfidenceMean(x,y,z) / dist(x,y,z);
            confidenceStd(x,y,z) = std::sqrt(std::max(march.getConfidenceStd(x,y,z) / dist(x,y,z) - confidenceMean(x,y,z)*confidenceMean(x,y,z),(T)0.0));
            confidenceMin(x,y,z) = march.getConfidenceMin(x,y,z);
        }
    }

    // Now post process the confidence minimum (save)
    cimg_mapXYZ(confidenceMin,x,y,z)
            confidenceMin(x,y,z) = (confidenceMin(x,y,z) == FLT_MAX) ? 0 : confidenceMin(x,y,z);

    if (cmean_name != "") {
        cout << "* Saving confidence mean to '" << cmean_name << "'... "; cout.flush();
        confidenceMean.save(cmean_name);
        cout << "OK" << endl;
    }

    if (cstd_name != "") {
        cout << "* Saving confidence std deviation to '" << cstd_name << "'... "; cout.flush();
        confidenceStd.save(cstd_name);
        cout << "OK" << endl;
    }

    if (cmin_name != "") {
        cout << "* Saving confidence minimum to '" << cmin_name << "'... "; cout.flush();
        confidenceMin.save(cmin_name);
        cout << "OK" << endl;
    }

    return 0;
}

