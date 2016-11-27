/* 
 * $Log: tools.c,v $
 * Revision 1.3  1995/11/09  12:59:28  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:20:12  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:20:12  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:58:19  ohenri
 * Initial revision
 * 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <malloc.h>
#include "tools.h"

/********************************************************/
int inside(int x,int y,int nx,int ny){return (x>=0 && y>=0 && x<nx && y<ny);}

float diffori(float x,float y)
{
  float xx,yy,diff;

  xx=mapangle(x,180);
  yy=mapangle(y,180);
  
  diff=MAX(xx,yy)-MIN(xx,yy);
  
  return ((diff > M_PI/2.0) ? (M_PI-diff):(diff));
}

float diffdir(float x,float y)
{
  float xx,yy,diff;

  xx=mapangle(x,360);
  yy=mapangle(y,360);
  
  diff=MAX(xx,yy)-MIN(xx,yy);
  
  return ((diff > M_PI) ? (2*M_PI-diff):(diff));
}



float mapangle(float x,int which)
{
  /* maps an angle x to [-pi,pi] or [0,pi] depending on which */
  
  while (x < -M_PI) x += 2*M_PI;
  while (x >  M_PI) x -= 2*M_PI;
  if ((which == 180)&&(x < 0)) x +=M_PI;

  return x; 
}

int angle2ori(float alpha,int nori)
{
  float incr,offset,angle;
  int i;

  angle=mapangle(alpha,180);
  incr=M_PI/nori;
  offset=incr/2.0;
  
  for (i=1;i<nori;i++){
    if ((angle <= (offset+i*incr))&&(angle > (offset+(i-1)*incr)))
      return i;
  }
  return 0;
}


int angle2discrete(float ang,int which)
{
  /* maps an arbitrary angle x to [0,1,2,3,4,5,6,7] or 
     [0,1,2,3] depending on which. */
  int no;
  float slice=M_PI/8.0;
  
  if ((ang > M_PI)||(ang < -M_PI)) ang=mapangle(ang,which);
  
  if ((ang > -slice)&&(ang < slice)) no=0;
  else if ((ang >= slice)&&(ang <= 3*slice))     no=1;
  else if ((ang > 3*slice)&&(ang < 5*slice))     no=2;
  else if ((ang >= 5*slice)&&(ang <= 7*slice))   no=3;
  else if ((ang >= -3*slice)&&(ang <= -slice))   no=7;
  else if ((ang > -5*slice)&&(ang < -3*slice))   no=6;
  else if ((ang >= -7*slice)&&(ang <= -5*slice)) no=5;
  else no=4;
  
  if ((which == 180)&&(no == 4)) no=0;

  return no;
}

/***********************************************/

float fbilin (float *buf, float x, float y, int nx)
{
  int ix, iy, pos;
  float rx, ry,tmp;

  ix = (int) x; iy = (int) y;
  rx = x - (float) ix; ry = y - (float) iy;
  tmp=rx*ry; pos=ix+iy*nx;

  return (buf [pos] * (1.0 - rx - ry + tmp) +
	  buf [pos+1] * (rx - tmp) +
	  buf [pos+nx] * (ry - tmp) +
	  buf [pos+nx+1] * tmp);
}

float cbilin (unsigned char *buf, float x, float y, int nx)
{
  int ix, iy, pos;
  float rx, ry,tmp;

  ix = (int) x; iy = (int) y;
  rx = x - (float) ix; ry = y - (float) iy;
  tmp=rx*ry; pos=ix+iy*nx;

  return ((float) buf [pos] * (1.0 - rx - ry + tmp) +
	  (float) buf [pos+1] * (rx - tmp) +
	  (float) buf [pos+nx] * (ry - tmp) +
	  (float) buf [pos+nx+1] * tmp);
}

/***********************************************/


int get_line(int x1,int y1,int x2,int y2,int **xret,int **yret)
{
  int d, x, y, ax, ay, sx, sy, dx, dy, ind, len;
  int dist=0;

  len=d4(x1,y1,x2,y2)+1;
  dx = x2-x1;  ax = ABS(dx)<<1;  sx = SGN(dx);
  dy = y2-y1;  ay = ABS(dy)<<1;  sy = SGN(dy);
  
  if (len > 1){
    *xret=(int *) malloc(len*sizeof(int));
    *yret=(int *) malloc(len*sizeof(int));
    
    x = x1;y = y1;
    if (ax>ay) {            
      d = ay-(ax>>1);
      for (;;) {
	(*xret)[dist]=x;(*yret)[dist++]=y;
	if (x==x2) return dist;
	if (d>=0) {y += sy;d -= ax;}
	x += sx;d += ay;
      }
    }
    else {
      d = ax-(ay>>1);
      for (;;) {
	(*xret)[dist]=x;(*yret)[dist++]=y;
	if (y==y2) return dist;
	if (d>=0) {x += sx;d -= ay;} 
	y += sy;d += ax;
      }
    }
  }
  else return 0;
}    

float de(int x1,int y1,int x2,int y2){return fsqrt(SQR(x1-x2)+SQR(y1-y2));}
int d8(int x1,int y1,int x2,int y2){return MAX(ABS(x1-x2),ABS(y1-y2));}
int d4(int x1,int y1,int x2,int y2){return ABS(x1-x2) + ABS(y1-y2);}

void draw_line(unsigned char *image,int nx,int ny,
	       int x1,int y1, int x2, int y2,unsigned char value)
{
  int *xx,*yy,dist,i;
  
  if (x1 < 0 || x1 >= nx || x2 < 0 || x2 >= nx ||
      y1 < 0 || y1 >= ny || y2 < 0 || y2 >= ny)
    fprintf(stderr,"draw_line : Illegal Coordinates\n");
  else{
    dist=get_line(x1,y1,x2,y2,&xx,&yy);
    for (i=0;i<dist;i++) image[xx[i]+yy[i]*nx]=value;
    free(xx);free(yy);
  }
}


/*****************************************************************************
 * int maxind (float *data, int n)
 *----------------------------------------------------------------------------
 * Function to get the index of the maximum value
 *****************************************************************************
 */

int maxind (float *data, int n)
{
  int i, maxptr = 0;
  float m = -1.0e38;

  for (i = 0; i < n; i++)
    if (data[i] > m) {m = data[i];maxptr = i;}
  
  return maxptr;
}

int minind (float *data, int n)
{
  int i, minptr = 0;
  float m;

  m = 1.0e38;
  for (i = 0; i < n; i++) {
    if (data[i] < m) {
      m = data[i];
      minptr = i;
    }
  }
  return (minptr);
}

static int mask3x3_x[] = {  1,  1,  1,  0,  0,  0, -1, -1, -1};
static int mask3x3_y[] = {  1,  0, -1,  1,  0, -1,  1,  0, -1};

int loc_max(float *map,int bd,int nx,int ny,
	    unsigned char *binary,float thresh)
{
  float level,array3x3[9];
  int x,y,yoff,i,pos,i_size=nx*ny,boarder=bd+1,mask[9],no=0;

  for (i=0;i<9;i++) mask[i]=mask3x3_x[i]+mask3x3_y[i]*nx; 

  memset ((void *) binary, 0, nx*ny*sizeof (char));
  level=find_absmax(map,bd,nx,ny)*thresh;
  
  for (y = boarder; y < (ny - boarder); y++) {
    yoff=y*nx;
    for (x = boarder; x < (nx - boarder); x++) {
      for (i = 0; i < 9; i++) 
	array3x3[i]=map[yoff+x+mask[i]];
      pos=maxind(array3x3,9);
      if ((pos == 4)&&(level<array3x3[4])){
	binary[yoff+x]=255;
	no++;
      }
    }
  }

  return no;
}
int loc_max_f(float *map,int bd,int nx,int ny,
	    float *fmap,float thresh)
{
  float level,array3x3[9];
  int x,y,yoff,i,pos,i_size=nx*ny,boarder=bd+1,mask[9],no=0;

  for (i=0;i<9;i++) mask[i]=mask3x3_x[i]+mask3x3_y[i]*nx; 

  memset ((void *) fmap, 0, nx*ny*sizeof (float));
  level=find_absmax(map,bd,nx,ny)*thresh;
  
  for (y = boarder; y < (ny - boarder); y++) {
    yoff=y*nx;
    for (x = boarder; x < (nx - boarder); x++) {
      for (i = 0; i < 9; i++) 
	array3x3[i]=map[yoff+x+mask[i]];
      pos=maxind(array3x3,9);
      if ((pos == 4)&&(level<array3x3[4])){
	fmap[yoff+x]=map[yoff+x];
	no++;
      }
    }
  }

  return no;
}

/*----------------------------------------------------------------------------
 * Smoothing a float image with a separable gaussian kernel.
 *****************************************************************************
 */

int smoothing (float *buf, int nx, int ny,float sigma)
{
  float *buf2,*mask,s;
  int support,i, j, k, kk, x, y;

  support = (int) 3.0*sigma;              /* set support size for gaussian */
  
  buf2=(float *) malloc(nx*ny*sizeof(float));
  mask=(float *) malloc((support+1)*sizeof(float));
  
  for(i=0; i <= support; i++)             /* one-sided gaussian mask */
    mask[i]=1.0/(sqrt(2.0*M_PI)*SQR(sigma))*fexp (-0.5*SQR((float)i / sigma));

  for(j = 0; j < nx*ny; j+=nx){           /* filter all lines */
    s = mask[0];
    for(i = j+support; i < j+nx-support; i++)     /* initialize vector with */
      buf2[i] = buf[i] * s;                       /* central weight */
    
    for(k=1; k<=support; k++){                   /* sum up remain. weighted */
      s = mask[k];                               /* vectors */
      for(i=j+support; i < j+nx-support; i++)                       
	buf2[i] += (buf[i+k] + buf[i-k]) * s;
    }
  }
  
  for(j = support*nx; j < (ny-support)*nx; j+=nx){ /* filter all columns */
    s = mask[0];
    for(i = j+support; i < j+nx-support; i++)      /* initialize vector */
      buf[i] = buf2[i] * s;                        /* central weight */
    
    for(k=1, kk=nx; kk<=support*nx; kk+=nx, k++){ /* sum up remaining  */
      s = mask[k];                                 /* weighted vectors */

      for(i=j+support; i < j+nx-support; i++)        
	buf[i] += (buf2[i+kk] + buf2[i-kk]) * s;
    }
  }

  for (i=0;i<support;i++){
    for (y=0;y<ny;y++){
      buf[i+y*nx]=0.0;
      buf[nx-i-1+y*nx]=0.0;
    }
    for (x=0;x<nx;x++){
      buf[x+i*nx]=0.0;
      buf[x+(ny-i-1)*nx]=0.0;
    }
  }
  free(mask);free(buf2);

  return support;
}



/*****************************************************************************
 * float find_absmax (float *buf, int nx, int ny)
 *----------------------------------------------------------------------------
 * Function to find max value in an float array.
 *****************************************************************************
 */
float find_absmax (float *buf,int boarder,int nx,int ny)
{
  int  yoff,x,y;
  float max=-1.0e38;
  
  for( y=boarder;y<(ny-boarder);y++){
    yoff=y*nx;
    for( x=boarder;x<(nx-boarder);x++){
      if (buf[x+yoff] > max) max=buf[x+yoff];
    }
  }

  return max;
}

int list_comp(Sortedlist *a,Sortedlist *b)
{
  Sortedlist aa=*a,bb=*b;
   
  if ((aa.value-bb.value) < 0.0) return 1;
  else return -1;
}

int generate_sortedlist(float *map,int nx,int ny,Sortedlist **list)
{
  int no=0,i,j,pos,i_size=nx*ny;
  
  for (i=0;i<i_size;i++) if (map[i] > 0.0) no++;

  *list=(Sortedlist *) malloc(no*sizeof(Sortedlist));
  
  for (j=no=0;j<ny;j++){
    for (i=0;i<nx;i++){
      pos=i+j*nx;
      if (map[pos] > 0.0){
	(*list)[no].x=i;
	(*list)[no].y=j;
	(*list)[no++].value=map[pos];
      }
    }
  }
  
  qsort((void *) *list,no,sizeof(Sortedlist),
	(int (*)(const void *,const void *)) list_comp);
  
  return no;
}

float Moments(float *xx,float *yy,int no)
{
  int i;
  float x_sum,y_sum,x2_sum,y2_sum,xy_sum,xm,ym;
  float mxx,mxy,myy;

  x_sum=y_sum=x2_sum=y2_sum=xy_sum=0.0;
  for (i=0;i<no;i++) {x_sum += xx[i];y_sum += yy[i];}
  xm = x_sum/(float) no;ym = y_sum/(float) no;
  mxx=myy=mxy=0.0;
  for (i=0;i<no;i++) {
    mxx += SQR(xx[i]-xm);
    myy += SQR(yy[i]-ym);
    mxy += (xx[i]-xm)*(yy[i]-ym);
  }

  return atan2((2.0 * mxy), (mxx - myy))/2.0;
}


int LSlinefit(float *x,float *y,int no,float r0_eps,float alpha0_eps,
	       float *a,float *b,float *c,float *error)
{
  int i,j,row,col,iterations=0;
  double alpha0,r0,cosa,sina,dalpha0,dr0,det;
  double N[2][2],L[2],Delta[2];
  double Ai[2],*w;
  
  /* get initial values */
  alpha0=mapangle(Moments(x,y,no)+M_PI*0.5,180);
  r0=x[0]*cos(alpha0)+y[0]*sin(alpha0);
  w=(double *) malloc(no*sizeof(double));
  
  do{
    cosa=cos(alpha0);sina=sin(alpha0);
    *error=0.0;
    for (i=0;i<no;i++) w[i]=x[i]*cosa+y[i]*sina-r0;
    
    L[0]=L[1]=N[0][0]=N[0][1]=N[1][0]=N[1][1]=0.0;
    for (i=0;i<no;i++){ 
      Ai[0]=-x[i]*sina+y[i]*cosa;Ai[1]=-1;
      N[0][0] += Ai[0]*Ai[0];
      N[0][1] += Ai[0]*Ai[1];
      N[1][0] += Ai[1]*Ai[0];
      N[1][1] += Ai[1]*Ai[1];
      L[0] += Ai[0]*w[i];
      L[1] += Ai[1]*w[i];
    }
    
    det=1.0/(N[0][0]*N[1][1]-N[0][1]*N[1][0]);
    dalpha0 = det*(N[1][1]*L[0] - N[0][1]*L[1]);
    dr0=det*(-N[1][0]*L[0] + N[0][0]*L[1]);
    alpha0 -= dalpha0;
    r0 -= dr0;
    iterations++;
  } while (fabs(dr0) > r0_eps || fabs(dalpha0) > alpha0_eps);
  
  
  free(w);
  *a=cos(alpha0);*b=sin(alpha0);*c= -r0;
  for (i=0,*error=0.0;i<no;i++) *error += fabs(x[i]*cosa+y[i]*sina-r0);
  *error /= (float) no;
  return iterations;
}

void ProjectPointOnLine(int x1,int y1,double a,double b,double c,
			float *xn1,float *yn1)
{
  double nx,ny,norm,x0,y0;

  if (fabs(b) > fabs(a)){x0=0;y0=-c/b;}
  else {y0=0;x0=-c/a;}

  norm=fsqrt(a*a+b*b);
  nx=-b/norm;ny=a/norm;
  *xn1 = x0 + ((double)(x1-x0)*nx+(double) (y1-y0)*ny)*nx;
  *yn1 = y0 + ((double)(x1-x0)*nx+(double) (y1-y0)*ny)*ny;
}

int InsideProjection(float xf,float yf,float *xf2,float *yf2)
{
  return (xf >= MIN(xf2[0],xf2[1]) && xf <= MAX(xf2[0],xf2[1]) &&
	  yf >= MIN(yf2[0],yf2[1]) && yf <= MAX(yf2[0],yf2[1]));
}


void sort(int n,float *ra)
{
  int l,j,ir,i;
  float rra;
  
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1)
      rra=ra[--l];
    else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
	ra[1]=rra;
	return;
      }
    }
    i=l;
    j=l << 1;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) ++j;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	j += (i=j);
      }
      else j=ir+1;
    }
    ra[i]=rra;
  }
}

/*----------------------------------------------------------------------------
 * Smoothing a 1D float array with a gaussian kernel.
 *****************************************************************************
 */

int smoothing1D (float *buf,int len,float sigma,float **buf2)
{
  float *mask,s;
  int support,i,k;
 
  support = (int) 3.0*sigma;                 /* set support size for gaussian */
  *buf2=(float *) malloc(len*sizeof(float));
  memset((void *) (*buf2),0,len*sizeof(float));
  mask=(float *) malloc((support+1)*sizeof(float));
    
  for(i=0; i <= support; i++)             /* one-sided gaussian mask */
    mask[i]=1.0/(sqrt(2.0*M_PI)*SQR(sigma))*fexp (-0.5*SQR((float)i / sigma));
  
  s = mask[0];
  for(i = support; i < len-support; i++)       /* initialize vector with */
    (*buf2)[i] = buf[i] * s;                      /* central weight */
  
  for(k=1; k<=support; k++){                   /* sum up remain. weighted */
    s = mask[k];                               /* vectors */
    for(i=support; i < len-support; i++)                       
      (*buf2)[i] += (buf[i+k] + buf[i-k]) * s;
  }
  
  free(mask);
  
  return support;
}

int make_histogram(float *buf,int nx,int bins,
		   float **hist,float *maxx,float *minx,float *incr)
{
  int i;

  *maxx=-1e38;
  *minx=1e38;
  
  for (i=0;i<nx;i++){
    if (buf[i] > *maxx) *maxx=buf[i];
    if (buf[i] < *minx) *minx=buf[i];
  }

  if (*maxx == *minx){*incr=0.0;return False;}
     
  *incr=(*maxx-*minx)/(float) bins;
  *hist=(float *) malloc(bins*sizeof(float));
  memset((void *) *hist,0,bins*sizeof(float));

  for (i=0;i<nx;i++) (*hist)[MIN(bins-1,(int) ((buf[i]-*minx)/ (*incr)))]++;
   
  return True;
}

void NRGB2HSV(float r,float g,float b,float *H,float *S,float *V)
{
  float x,v,s,R1,G1,B1;

  /* convert NRGB to HSV */

  v=*V=MAX(r,MAX(g,b));
  x=MIN(r,MIN(g,b));
  s=*S=(v - x)/v;
  if (s == 0.0) *H=0.0;
  else{
    R1=(v-r)/(v-x);
    G1=(v-g)/(v-x);
    B1=(v-b)/(v-x);
    if (r == v){
      if (g == x) *H=5.0+B1;
      else *H=1.0-G1;
    }
    else if (g == v){
      if (b == x) *H=R1+1.0;
      else *H=3.0-B1;
    }
    else {
      if (r == x) *H=3.0+G1;
      else *H=5.0-R1;
    }
  }
   
  (*H) *= 60.0*GRa; /* transfer hexagonal coordinates to [0.0 -> 2*M_PI]  */
  (*S) *= 100.0;    /* transfer saturation [0.0 -> 1.0] to [0.0 -> 100.0] */
  (*V) *= 100.0;    /* transfer value [0.0 -> 1.0] to [0.0 -> 100.0]      */
}

void HSV2NRGB(float H,float S,float V,float *r,float *g,float *b)
{
  float H1,F,A[7];
  int I;

  H1 = H/(60.0*GRa); /* transfer hue [0.0 -> 2*M_PI] to hexagonal coordinates */
  S /= 100.0;        /* transfer saturation  [0.0 -> 100.0] to [0.0 -> 1.0]   */
  V /= 100.0;        /* transfer value [0.0 -> 100.0] to [0.0 -> 1.0]         */
  I=(int) H1;
  F=H1-I;
  A[1]=A[2]=V;
  A[3]=V*(1.0-(S*F));
  A[4]=A[5]=V*(1.0-S);
  A[6]=V*(1.0-(S*(1.0-F)));
  if (I>4) I -= 4;
  else I += 2;
  *r=A[I];
  if (I>4)  I -= 4;
  else I += 2;
  *b=A[I];
  if (I>4)  I -= 4;
  else I += 2;
  *g=A[I];
}

void RGB2Lab(float R,float G,float B,float *L,float *a,float *b)
{
  float X,Y,Z,X0,Y0,Z0;
  float smallConst=16.0/116.0;
  float quotientThresh=0.008856;

  /* Celenk (Computer Vision, Graphics and Image processing 52, 145-170 ) */
  X= 2.7690*R + 1.7518*G + 1.1300*B;
  Y= 1.0000*R + 4.5907*G + 0.0601*B;
  Z= 0.0000*R + 0.0565*G + 5.5943*B;
  
  /* float Rwhite=255.0,Gwhite=255.0,Bwhite=255.0;
     X0=2.7690*Rwhite + 1.7518*Gwhite + 1.1300*Bwhite; = 5.6508 * 255 = 1442.7900
     Y0=1.0000*Rwhite + 4.5907*Gwhite + 0.0601*Bwhite; = 5.6508 * 255 = 1442.7900
     Z0=0.0000*Rwhite + 0.0565*Gwhite + 5.5943*Bwhite; = 5.6508 * 255 = 1442.7900
     
     With RGB values between 0 .. 255 we achieve the following [min,max] values:

     Lmin=0.000000     Lmax=99.950775
     amin=-128.091721  amax=182.353058
     bmin=-155.282272  bmax=156.135483
     */
  
  X /= 1442.7900;
  Y /= 1442.7900;
  Z /= 1442.7900;
  
  /*  preparation of transform (Wyszecki and Stiles p. 167) */
  if (Y <= quotientThresh) {
    *L = 903.3 * Y;
    Y = 7.787 * Y + smallConst;
  }
  else {
    Y = fpow(Y,1.0/3.0);
    *L= 116.0*Y - 16.0;
  }
  
  if (X > quotientThresh) X = fpow(X,1.0/3.0);
  else X = 7.787 * X + smallConst;
  
  if (Z > quotientThresh) Z = fpow(Z,1.0/3.0);
  else Z = 7.787 * Z + smallConst;
  
  *a= 500*(X - Y);
  *b= 200*(Y - Z);
}


void RGB2XYZ(float R,float G,float B,float *X,float *Y,float *Z)
{
  *X= 2.7690*R + 1.7518*G + 1.1300*B;
  *Y= 1.0000*R + 4.5907*G + 0.0601*B;
  *Z= 0.0000*R + 0.0565*G + 5.5943*B;
}


void RGB2Luv(float R,float G,float B,float *L,float *u,float *v)
{
  float X,Y,Z,Xn,Yn,Zn,up,vp,upn,vpn;
  float smallConst=16.0/116.0;
  float quotientThresh=0.008856;

  /* Celenk (Computer Vision, Graphics and Image processing 52, 145-170 ) */
  X= 2.7690*R + 1.7518*G + 1.1300*B;
  Y= 1.0000*R + 4.5907*G + 0.0601*B;
  Z= 0.0000*R + 0.0565*G + 5.5943*B;
  
  Xn = 1442.7900;
  Yn = 1442.7900;
  Zn = 1442.7900;
  
  if (Y/Yn <= quotientThresh) *L = 903.3 * Y/Yn;
  else *L= 116.0*fpow(Y/Yn,1.0/3.0) - 16.0;
 
  up=4.0*X / (X + 15.0*Y + 3.0*Z);
  vp=9.0*Y / (X + 15.0*Y + 3.0*Z);
  upn=4.0*Xn / (Xn + 15.0*Yn + 3.0*Zn);
  vpn=9.0*Yn / (Xn + 15.0*Yn + 3.0*Zn);
     
  *u= 13*(*L)*(up-upn);
  *v= 13*(*L)*(vp-vpn);
}



/****************************************************************************/
void sofit (float *buf, int x, int y, int nx, float *a)

/*--------------------------------------------------------------------------*/
/* 2D second order fit in 3x3 neighborhood                                  */
/* f(x,y) = a0 + a1 x * a2 y + a3 x^2 + a4 y^2 + a5 xy                      */
/* 3x3 surround index definition and offsets                                */
/*     ------------                                                         */
/*    | 4 | 3 | 2 |    x-1+nx*(y+1)  |  x  +nx*(y+1)  |  x+1+nx*(y+1)       */
/*    | 5 | 0 | 1 |    x-1+nx*y      |  x  +nx*y      |  x+1+nx*y           */
/*    | 6 | 7 | 8 |    x-1+nx*(y-1)  |  x  +nx*(y-1)  |  x+1+nx*(y-1)       */
/*    ------------                                                          */

{
  float sxp, sxn, syp, syn, stot;

  sxp = buf[x+1+nx* y   ] + buf[x+1+nx*(y+1)] + buf[x+1+nx*(y-1)]; /* 1,2,8 */
  sxn = buf[x-1+nx*(y+1)] + buf[x-1+nx* y   ] + buf[x-1+nx*(y-1)]; /* 4,5,6 */
  syp = buf[x+1+nx*(y+1)] + buf[x  +nx*(y+1)] + buf[x-1+nx*(y+1)]; /* 2,3,4 */
  syn = buf[x-1+nx*(y-1)] + buf[x  +nx*(y-1)] + buf[x+1+nx*(y-1)]; /* 6,7,8 */
  stot  =  buf[x  +nx*y] + buf[x  +nx*(y+1)] + buf[x  +nx*(y-1)]   /* 0,3,7 */
         + sxp + sxn;

  a[0] = 0.5555555555 * stot - (sxn+sxp+syn+syp)/3.0; 
  a[1] = (sxp - sxn)/6.0;
  a[2] = (syp - syn)/6.0;
  a[3] = (sxp + sxn)/2.0 - stot/3.0;
  a[4] = (syp + syn)/2.0 - stot/3.0;
  a[5] = (  buf[x+1+nx*(y+1)] + buf[x-1+nx*(y-1)]        /* 2,6 */
          - buf[x-1+nx*(y+1)] - buf[x+1+nx*(y-1)])/4.0;  /* 4,8 */

  return;
}


/****************************************************************************/
void sofit2 (float *buf, int x, int y, int nx, float *a)

/*--------------------------------------------------------------------------*/
/* 2D second order fit in 5x5 neighborhood                                  */
/* f(x,y) = a0 + a1 x * a2 y + a3 x^2 + a4 y^2 + a5 xy                      */
/* 5x5 surround index definition and offsets                                */
/*  -------------------------                                               */
/* | 15 | 14 | 13 | 12 | 11 |                                               */
/* | 16 |  4 |  3 |  2 | 10 |                                               */ 
/* | 17 |  5 |  0 |  1 |  9 |          Matrix Indices                       */
/* | 18 |  6 |  7 |  8 | 24 |                                               */
/* | 19 | 20 | 21 | 22 | 23 |                                               */
/* -------------------------                                                */
/* Constants:                                                               */
/* 27/175 = 0.1542857143                                                    */
/* 1/35   = 0.0285714286                                                    */
/* 1/50   = 0.02                                                            */
/* 1/70   = 0.0142857143                                                    */

{
  float sx0, sxp1, sxn1, syp1, syn1,  sxp2, sxn2, syp2, syn2, stot;



  sx0  = buf[x  +nx* y   ]+buf[x  +nx*(y+1)]+buf[x  +nx*(y-1)] + /* 0,3,7   */
         buf[x  +nx*(y+2)]+buf[x  +nx*(y-2)];                    /* 13,21   */

  sxp1 = buf[x+1+nx* y   ]+buf[x+1+nx*(y+1)]+buf[x+1+nx*(y-1)] + /* 1,2,8   */
         buf[x+1+nx*(y+2)]+buf[x+1+nx*(y-2)];                    /* 12,22   */
  sxn1 = buf[x-1+nx*(y+1)]+buf[x-1+nx* y   ]+buf[x-1+nx*(y-1)] + /* 4,5,6   */
         buf[x-1+nx*(y+2)]+buf[x-1+nx*(y-2)];                    /* 14,20   */
  sxp2 = buf[x+2+nx*(y+2)]+buf[x+2+nx*(y+1)]+buf[x+2+nx* y   ] + /* 11,10,9 */
         buf[x+2+nx*(y-1)]+buf[x+2+nx*(y-2)];                    /* 24,23   */
  sxn2 = buf[x-2+nx*(y+2)]+buf[x-2+nx*(y+1)]+buf[x-2+nx* y   ] + /* 15,16,17*/
         buf[x-2+nx*(y-1)]+buf[x-2+nx*(y-2)];                    /* 18,19   */

  syp1 = buf[x-1+nx*(y+1)]+buf[x  +nx*(y+1)]+buf[x+1+nx*(y+1)] + /* 4,3,2   */
         buf[x-2+nx*(y+1)]+buf[x+2+nx*(y+1)];                    /* 16,10   */
  syn1 = buf[x-1+nx*(y-1)]+buf[x  +nx*(y-1)]+buf[x+1+nx*(y-1)] + /* 6,7,8   */
         buf[x-2+nx*(y-1)]+buf[x+2+nx*(y-1)];                    /* 18,24   */
  syp2 = buf[x-2+nx*(y+2)]+buf[x-1+nx*(y+2)]+buf[x  +nx*(y+2)] + /* 15,14,13*/
         buf[x+1+nx*(y+2)]+buf[x+2+nx*(y+2)];                    /* 12,11   */
  syn2 = buf[x-2+nx*(y-2)]+buf[x-1+nx*(y-2)]+buf[x  +nx*(y-2)] + /* 19,20,21*/
         buf[x+1+nx*(y-2)]+buf[x+2+nx*(y-2)];                    /* 22,23   */

  stot = sx0 + sxp1 + sxp2 + sxn1 + sxn2;

  a[0] =  0.1542857143 * stot - 
          0.0285714286 * (sxn1+sxp1+syn1+syp1 + 4.0*(sxn2+sxp2+syn2+syp2)); 
  a[1] =  0.02 * (sxp1 - sxn1 + 2.0*(sxp2 - sxn2));
  a[2] =  0.02 * (syp1 - syn1 + 2.0*(syp2 - syn2));
  a[3] =  0.0142857143 * (sxp1 + sxn1 + 4.0*(sxp2 + sxn2)) - 
          0.0285714286 * stot;
  a[4] =  0.0142857143 * (syp1 + syn1 + 4.0*(syp2 + syn2)) - 
          0.0285714286 * stot;
  a[5] = 0.01 * (4.0 *(buf[x-2+nx*(y-2)] + buf[x+2+nx*(y+2)] -  /* 19,11 */
                       buf[x-2+nx*(y+2)] - buf[x+2+nx*(y-2)]) + /* 15,23 */
                 2.0 *(buf[x-2+nx*(y-1)] + buf[x-1+nx*(y-2)] +  /* 18,20 */
                       buf[x+1+nx*(y+2)] + buf[x+2+nx*(y+1)] -  /* 12,10 */
		       buf[x-1+nx*(y+2)] - buf[x-2+nx*(y+1)] -  /* 14,16 */
		       buf[x+1+nx*(y-2)] - buf[x+2+nx*(y-1)]) + /* 22,24 */
		       buf[x+1+nx*(y+1)] + buf[x-1+nx*(y-1)] -  /*  2,6 */
		       buf[x-1+nx*(y+1)] - buf[x+1+nx*(y-1)]);  /*  4,8 */

  return;
}
/*==========================================================================*/


