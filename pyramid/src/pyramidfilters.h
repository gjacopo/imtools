/* ----------------------------------------------------------------------------
	Filename:  	pyramidfilters.h
	
	Project:	Biomedical Imaging Library
 
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		17 March 1999
	
	Purpose:	Header associated to pyramidfilters.c
			 
---------------------------------------------------------------------------- */
 void PyramidFilterSplinel2(double g[],long *ng,double *h,long *nh,long Order);
 void PyramidFilterSplineL2(double g[],long *ng,double *h,long *nh,long Order);
 void PyramidFilterCentered(double g[],long *ng,double h[],long *nh,long Order);
 void PyramidFilterCenteredL2(double g[],long *ng,double h[],long *nh,long Order);
 void PyramidFilterCenteredL2Derivate(double g[],long *ng,double h[],long *nh,long Order);

