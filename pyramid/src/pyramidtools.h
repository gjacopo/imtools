/* ----------------------------------------------------------------------------
	Filename:  	pyramidtools.h
	
	Project:	Biomedical Imaging Library
 double
	Author:		Daniel Sage
				Swiss Federal Institute of Technology - Lausanne
				Biomedical Imaging Group
				EPFL/DMT/IOA, BM-Ecublens, CH-1015 Lausanne, Switzerland
				
	Date:		17 March 1999
	
	Purpose:	Header associated to pyramidtools.c
			 
---------------------------------------------------------------------------- */

 int GetPyramidFilter(
				char *Filter, 
				long Order, 
				double g[],long *ng,
				double h[],long *nh, 
				short *FlagCentered		
				);	
					
 int Reduce_2D(	
				double *In, long NxIn, long NyIn,
				double *Out,
				double w[], long nw,
				short FlagCentered
				);
				
 int Expand_2D(	
				double *In, long NxIn, long NyIn,
				double *Out,
				double w[], long nw,
				short FlagCentered
				);

 void Reduce_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
 void Expand_1D(	
				double x[], long nx,
				double y[],
				double w[], long nw,
				short FlagCentered
				);
						
