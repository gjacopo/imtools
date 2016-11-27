/* 
 * $Log: tools.h,v $
 * Revision 1.2  1995/10/02  12:20:52  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:20:52  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:50:13  ohenri
 * Initial revision
 * 
 *
 */

#if !defined(_OHENRI_TOOLS_)
#define _OHENRI_TOOLS_

#if !defined(_OHENRI_DEFS_)
#include "defs.h"
#endif

#define NCMAP 11
#define RAS_MAGIC       0x59a66a95
#define RT_STANDARD     1       /* Raw pixrect image in 68000 byte order */

                                /* Sun supported ras_maptype's */
#define RMT_NONE        0       /* ras_maplength is expected to be 0 */
#define RMT_EQUAL_RGB   1       /* red[ras_maplength/3],green[],blue[] */
#define FLOAT 4
#define CHAR  1

typedef struct {int x,y;float value;} Sortedlist;

typedef int Koord;
typedef struct {Koord x, y;} Punkt;
typedef struct {Punkt von, nach;} Linie;
typedef float Trafo[2][3];
typedef struct {int length; Punkt *pkt;} Polygon; 
typedef struct {int count; Polygon *poly;} Polygon_list;
typedef struct {int clsize;unsigned int clid;int rank;} cldsc;

/*********************************************************
 *	pgon_approx.c
 *********************************************************/
int polygon_approx(Punkt *,int,float,float,Punkt *);
int polygon_approx_i(Punkt *,int,float,float,Punkt *);
float Tolerance(float,float);
/* void SetTolerance(float,float,float,float); */

/*********************************************************
 *	thin.c
 *********************************************************/
int thin(unsigned char *,int,int,unsigned char *);
static int friss(int *);
static void thinlut(unsigned char *);
int elemcount(int);
static init(unsigned char *,int,int,unsigned char *);

/*********************************************************
 *      ccl.c
 *********************************************************/
int ccl(unsigned char *,unsigned int *,int,int,int,unsigned int **);
int clcomp(cldsc *,cldsc *);
int clback(cldsc *,cldsc *);
int clmerge(unsigned int *,int *,int,int,int);
int clanal(int,unsigned char *,unsigned int *,int,int,unsigned int *);

/*********************************************************
 *	tools.c
 *********************************************************/
int inside(int,int,int,int);
float diffori(float,float);
float diffdir(float,float);
float mapangle(float,int);
int angle2ori(float,int);
int angle2discrete(float,int);
float fbilin (float *, float, float, int);
float cbilin (unsigned char *, float, float, int);
int get_line(int,int,int,int,int **,int **);
float de(int,int,int,int);
int d8(int,int,int,int);
int d4(int,int,int,int);
void draw_line(unsigned char *,int,int,int,int,int,int,unsigned char);

int maxind (float *, int);
int minind (float *, int);
float find_absmax (float *,int,int,int);
int smoothing (float *, int, int,float);
int loc_max(float *,int,int,int,unsigned char *,float);

int list_comp(Sortedlist *,Sortedlist *);
int generate_sortedlist(float *,int ,int ,Sortedlist **);

float Moments(float *,float *,int);
int LSlinefit(float *,float *,int,float,float,float *,float *,float *,float *);
void ProjectPointOnLine(int,int,double,double,double,float *,float *);
int InsideProjection(float,float,float *,float *);
void sort(int,float *);
int smoothing1D (float *,int,float,float **);
int make_histogram(float *,int,int,float **,float *,float *,float *);
void NRGB2HSV(float,float,float,float *,float *,float *);
void HSV2NRGB(float,float,float,float *,float *,float *);
void RGB2Lab(float,float,float,float *,float *,float *);
void RGB2XYZ(float,float,float,float *,float *,float *);
void RGB2Luv(float,float,float,float *,float *,float *);
void sofit (float *,int,int,int,float *);
void sofit2 (float *,int,int,int,float *);

/*********************************************************
 *	io.c
 *********************************************************/
int read_image(unsigned char **,char *,int *,int *);
int read_RGB_image(unsigned char **,unsigned char **,unsigned char **,
		   char *,int *,int *);
int write_grey_image(void *,char *,int,int,int,int);
int write_RGB_image(unsigned char *,unsigned char *,unsigned char *,
		    char *,int,int);
void write_color_image(unsigned char *,char *,int,int);
void write_overlay_image(unsigned char *,unsigned char *,char *,int,int,
			 char,char,char);

char *get_generic(char *);
float get_float_argument(int,char **,char *);
char *get_char_argument(int,char **,char *);
int exists_argument(int,char **,char *);

/*********************************************************
 *	useful.c
 *********************************************************/
void normalize(float *,int,int);
void float_nms(float *,int,int,int,unsigned char *,float,float *);
void maxvalues(float *,int,int,unsigned char *,float *,float *,float *);

/*********************************************************
 *	fill.c
 *********************************************************/
void Fill(int,int,unsigned char *,int,int,unsigned charr);



#endif
