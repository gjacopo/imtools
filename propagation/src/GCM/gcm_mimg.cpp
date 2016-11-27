/*------------------------------------------------------------------------------------------------------
  
  File        : gcm_MImg.cpp

  Description : Example of use of the GCM library

  Author      : Emmanuel Prados (UCLA/INRIA), Christophe Lenglet (INRIA), Jean-Philippe Pons (INRIA)
  
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
  
  --------------

  Associated publications  : This C++ code corresponds to the implementation of the algorithm presented in the following articles:  
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; Control Theory and Fast Marching Methods for Brain Connectivity Mapping; INRIA Research Report 5845 -- UCLA Computer Science Department Technical Report 060004, February 2006.
  - E. Prados, C. Lenglet, J.P. Pons, N. Wotawa, R. Deriche, O. Faugeras, S. Soatto; Control Theory and Fast Marching Methods for Brain Connectivity Mapping; Proc. IEEE Computer Society Conference on Computer Vision and Pattern Recognition, New York, NY, I: 1076-1083, June 17-22, 2006.  

  Please, if you use the GCM library in you work, make sure you will include the reference to the work of the authors in your publications.
  
  --------------
  
  Technical Comments and Detailled description:
   - Goal: example of use of the GCM library.
   - Tasks:
   * loading of the data (image, starting point, optional mask)
   * initialization of the objects
   * Computation the distance function/optimal dynamics/confidence statistics within the mask
   * Saving distance function/optimal dynamics/confidence statistics
  
  ----------------------------------------------------------------------------------------------------*/



#include <iostream>
#include <cfloat>
#include <climits>
#include <vector>
#include <string>
#include "MImg.h"
#include "AnisotropicTensorDistance.h"
#include "AnisotropicTensorDistanceConfidence.h"
using namespace cimg_library;
using namespace std;
using namespace FastLevelSet;

typedef float T;

template<typename T>
std::vector<MImg<T> > singlePointConnectivity(const MImg<T>& tensor, const MImg<T>& tensor_power, const MImg<T>& mask, const int& x0, const int& y0, const int& z0, const T& alpha) {
    std::vector<MImg<T> > results; // Distance, dynamics, mean, std, min

    // Fast marching initialization
    cout << std::endl; cout.flush();
    cout << "* Initializing for voxel " << x0 << ", " << y0 << ", " << z0 << "... "; cout.flush();
    MImg<T> dist(tensor.width,tensor.height,tensor.depth);
    dist.fill(FLT_MAX);
    AnisotropicTensorDistanceConfidence<T> march(dist.data,dist.width,dist.height,dist.depth,mask.data,tensor.data,tensor_power.data,alpha);
    cout << "OK" << endl;

    dist(x0,y0,z0) = 0;
    march.AddTrialPoint(x0,y0,z0);

    // Compute the distance function, the optimal dynamics and confidence statistics, within the mask
    cout << "* Working... "; cout.flush();
    march.Run();
    cout << "OK" << endl;

    // Store it
    results.push_back(dist);

    // Extract optimal dynamics
    // It is already computed; the computational result of the optimal dynamics is stacked in some variables of the "march" object)
    MImg<T> optDynamics(tensor.width,tensor.height,tensor.depth,3);
    cimg_mapXYZ(optDynamics,x,y,z) {
        optDynamics(x,y,z,0) = march.getOptDynamics(x,y,z,0);
        optDynamics(x,y,z,1) = march.getOptDynamics(x,y,z,1);
        optDynamics(x,y,z,2) = march.getOptDynamics(x,y,z,2);
    }

    // Store it
    results.push_back(optDynamics);

    // Get confidence
    // As the optimal dynamics, the confidence is already computed)
    MImg<T> confidenceMean(tensor.width,tensor.height,tensor.depth);
    MImg<T> confidenceStd(tensor.width,tensor.height,tensor.depth);
    MImg<T> confidenceMin(tensor.width,tensor.height,tensor.depth);

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

    // Now post process the confidence minimum
    cimg_mapXYZ(confidenceMin,x,y,z)
            confidenceMin(x,y,z) = (confidenceMin(x,y,z) == FLT_MAX) ? 0 : confidenceMin(x,y,z);

    // Store everything
    results.push_back(confidenceMean);
    results.push_back(confidenceStd);
    results.push_back(confidenceMin);

    return results;
}


int main(int argc, char **argv) {

    // Parameters
    cimg_usage("Geodesic Connectivity Mapping (Index of connectivity)\n");

    const string tensors_name = cimg_option("-i","","Diffusion Tensor image : filein DTI");
    const string origin = cimg_option("-p","-1,-1,-1","Point of interest : string");
    const string roi = cimg_option("-r","","Region of interest : filein mask | OPTIONAL");
    const string mask_name = cimg_option("-m","","Mask : filein mask | OPTIONAL");
    const T alpha = cimg_option("-a",0.0,"Metric power : float");
    const string distance_name = cimg_option("-dist","distance.hdr","Distance output : fileout SCL");
    const string odyn_name = cimg_option("-odyn","dynamics.hdr","Optimal dynamics output : fileout SCL");
    const string cmean_name = cimg_option("-mean","connectivity_mean.hdr","Confidence mean output : fileout SCL");
    const string cstd_name = cimg_option("-std","connectivity_std.hdr","Confidence std output : fileout SCL");
    const string cmin_name = cimg_option("-min","connectivity_min.hdr","Confidence min output : fileout SCL");

    if (cimg_option("-h",(const char *)0,0)) return 0;
    if (argc<4) {
        cerr << "Not enough arguments" << endl;
        return 1;
    }

    // Load the tensor field
    cout << "* Loading tensor field from '" << tensors_name << "'... "; cout.flush();
    if (tensors_name == "") throw CImgException("GCM: Please provide a Diffusion Tensor Image");
    MImg<T> tensor(tensors_name.c_str());
    cout << tensor.width << " x " << tensor.height << " x " << tensor.depth << " x " << tensor.dim << " "; cout.flush();
    cout << "OK" << endl;

    // Load the mask
    cout << "* Loading mask from '" << mask_name << "'... "; cout.flush();
    if (mask_name == "") throw CImgException("GCM: Please provide a mask of the region (full brain, hemisphere, segmentation ...) to work on");
    MImg<T> mask(mask_name.c_str());
    cout << mask.width << " x " << mask.height << " x " << mask.depth << " "; cout.flush();
    cout << "OK" << endl;

    // Tensors power
    MImg<T> tensor_power = tensor;

    // tensor_power already contains tensor and if alpha==0 the code in
    // AnisotropicTensorDistanceConfidence.h (l. 62) will explicitely compute the Euclidean norm of the dynamics
    if (alpha != 0.0 && alpha != 1.0) {
        cout << "* Computing tensors power for confidence estimation ... "; cout.flush();

        MImg<T> val(3,1,1,1), valM(3,3,1,1), vec(3,3,1,1);
        valM.fill(0.0);
		vec.fill(0.0);

        cimg_mapXYZ(tensor_power,x,y,z) {
            const MImg<T> D = tensor_power.get_tensor(x,y,z);
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

    // Compute everything
    MImg<T> dummyinit(tensor.width, tensor.height, tensor.depth);
    dummyinit.fill(0);
    std::vector<MImg<T> > results;

    int x0,y0,z0;
    if (origin != "-1,-1,-1") {
        sscanf(origin.c_str(),"%d,%d,%d",&x0,&y0,&z0);
        results = singlePointConnectivity(tensor, tensor_power, mask, x0, y0, z0, alpha);
    } else if (roi != "") {
        MImg<T> seeds(roi.c_str());
        for (unsigned i=0; i<5; ++i) results.push_back(dummyinit);

        cimg_mapXYZ(seeds,x,y,z) {
            if (seeds(x,y,z)>0) {
                std::vector<MImg<T> > voxelresults = singlePointConnectivity(tensor, tensor_power, mask, x, y, z, alpha);
                for (unsigned i=0; i<5; ++i) results[i] += voxelresults[i];
            }
        }
    } else throw CImgException("GCM: Please provide a point or a region of interest");

    cout << std::endl; cout.flush();

    // Save the result
    if (distance_name != "") {
        cout << "* Saving distance to '" << distance_name << "'... "; cout.flush();
        results[0].SetTransform(tensor.GetTransform());
        results[0].Save(distance_name.c_str());
        cout << "OK" << endl;
    }

    if (odyn_name != "") {
        cout << "* Saving optimal dynamics to '" << odyn_name << "'... "; cout.flush();
        results[1].SetTransform(tensor.GetTransform());
        results[1].Save(odyn_name.c_str());
        cout << "OK" << endl;
    }

    if (cmean_name != "") {
        cout << "* Saving connectivity mean to '" << cmean_name << "'... "; cout.flush();
        results[2].SetTransform(tensor.GetTransform());
        results[2].Save(cmean_name.c_str());
        cout << "OK" << endl;
    }

    if (cstd_name != "") {
        cout << "* Saving connectivity standard deviation to '" << cstd_name << "'... "; cout.flush();
        results[3].SetTransform(tensor.GetTransform());
        results[3].Save(cstd_name.c_str());
        cout << "OK" << endl;
    }

    if (cmin_name != "") {
        cout << "* Saving confidence minimum to '" << cmin_name << "'... "; cout.flush();
        results[4].SetTransform(tensor.GetTransform());
        results[4].Save(cmin_name.c_str());
        cout << "OK" << endl;
    }

    return 0;
}

