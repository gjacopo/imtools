/* 
 * $Log: lowlevel.h,v $
 * Revision 1.2  1995/10/02  12:23:59  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:50:13  ohenri
 * Initial revision
 *
 * Revision 1.1  1995/06/21  12:50:13  ohenri
 * Initial revision
 * 
 *
 */

#if !defined(_OHENRI_LOWLEVEL_)
#define _OHENRI_LOWLEVEL_

#if !defined(_OHENRI_TOOLS_)
#include "tools.h"
#endif

#if !defined(_OHENRI_DEFS_)
#include "defs.h"
#endif

/*********************************************************
 *	egabor.c
 *********************************************************/
void egabor_filter(float *,int,int,float,float,float,int,float,float *,
		   float *,float *);
float egabor_params(float,float,float,float,int,int);
int embed(float *,int,int,float **,int *,int *);
void crop(float *,int,int,int,int);
float *cvrtFloatToFourier (float *,int,int);
float *cvrtFourierToFloat (float *,int,int);
void four1 (float *,int,int);void realft (float *,int ,int);
void fft_alloc (int,int);
void fft_free (void);
void fft (float *,float *,int);

/*********************************************************
 *     noise_estimation.c
 *********************************************************/
void NoiseEstimation(float *,int,int,float,int,float *,float *,float *,float *);

/*********************************************************
 *     contour_detection.c
 *********************************************************/
#define NULLTYPE     0
#define EDGETYPE   200
#define LINETYPE   100

void contour_detection(float *,float *,float *,float,int,int,int,float *,
		       float *,float *,unsigned char *,float,float,float);
void sofit (float *, int , int , int, float *);
void sofit2 (float *, int , int , int, float *);


/*********************************************************
 *	kpt.c
 *********************************************************/
int keypoint_map(float *,int,int,int,float,float,float,float,float *);

/*********************************************************
 *	maps.c
 *********************************************************/
void quality_dft(float *,float *, int,int,float,float *);
void amplitude_phase(float *,int,int,float *,float *,int);

/*********************************************************
 *	start.c
 *********************************************************/
void Start_points(unsigned char *,float *,float *,float *,float *,
		  unsigned char *,float *,int,int,int,float);

/*********************************************************
 *	ll_processing.c
 *********************************************************/
void Low_level_Processing(unsigned char *,int,int,int,
			  float,float,float,float,float,float,
			  float,float,float,float,
			  float *,float *,float *,
			  float *,unsigned char *,
			  float *,float *,int);

/*********************************************************
 *	canny.c
 *********************************************************/
int canny(unsigned char *,float *,float *,float *,int,int,float);
     
#endif
