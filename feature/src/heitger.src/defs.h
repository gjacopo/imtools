/*
 * $Log: defs.h,v $
 * Revision 1.4  1995/10/02  12:23:59  ohenri
 * *** empty log message ***
 *
 * Revision 1.3  1995/06/21  14:48:31  ohenri
 * *** empty log message ***
 *
 * Revision 1.3  1995/06/21  14:48:31  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/06/21  13:25:44  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:35:01  ohenri
 * Initial revision
 *
 * Revision 1.1  1995/06/21  12:35:01  ohenri
 * Initial revision
 *
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>

#if !defined(_OHENRI_DEFS_)
#define _OHENRI_DEFS_

#if !defined( __sgi)
#define fexp(x) ((float) exp ((double) x))
#define flog(x) ((float) log10 ((double)(x)))
#define fcos(x) ((float) cos ((double) x))
#define fsin(x) ((float) sin ((double) x))
#define ftan(x) ((float) tan ((double) x))
#define fsqrt(x) ((float) sqrt ((double) x))
#define fatan2(x,y) ((float) atan2((double)(x),(double)(y)))
#define fpow(x,y) ((float) pow ((double) (x), (double) (y)))
#endif

#define ffabs(x)  ((float) fabs ((double) (x)))
#define SQR(x) ((x)*(x))
#define ABS(x) (((x) > 0) ? (x) : (-1*(x)))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define SGN(a)   (((a)<0) ? -1 : 1)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* FBI = linear interpolation */
#define FBI(buf,d) ((buf)[(int)(d)]+((buf)[(int)(d)+1]-(buf)[(int)(d)])*((d)-(float)(int)(d)))

#define FALSE 0
#define TRUE 1
#define True 1
#define False 0
#define IROUND(x) (((x) > 0) ? (int)((x) + 0.5) : (int)((x) - 0.5))
#define CLIP(x) (((x) > 0) ? (x) : 0.0)

#define PPi    57.2957795
#define GRa    0.0174532925
#define LARGE 1e10
#define VERY_LARGE 1e38
#define VERY_SMALL 1e-38

#endif
