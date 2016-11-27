/*===============================================================
 *
 * Institute of Communications Technology, Image Vision Group
 * Swiss Federal Institute of Technology at Zurich, Switzerland.
 *
 *===============================================================
 */


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <memory.h>
#include <malloc.h>
#include "lowlevel.h"

static char errmessage[] = "Usage:  Keypoints -i generic_filename [options]\n        ,where options are any of the following\n\n        -thresh [float]   : threshold (in %%) relative maximum (def. 20.0)\n        -size [int]       : support for local maximum, e.g. 3,5,7,... (def. 3)\n\n        -ascii            : output is an ascii file with [x,y,magnitude]\n        -overlay [mode]   : output is an (color) overlaid raster file\n                            where mode is either cross or circle\n        -back filename    : background image for overlay (rasterfile)\n        -radius [float]   : size of overlay marker (def. 2.5)\n\n        The program performs a local maximum operation on the keypoint map\n        and a thresholding operation. The generic filename must be the same\n        as the one used after -o for the SE program!\n";

void main(int argc,char **argv)
{
  char *fgen=NULL,fout[256],*fovl,arr[256],*type;
  int i,d,nx,ny,s,size=3,nms,x,y,dx,dy,nx2,ny2;
  float thresh=20.0;
  int ascii=False,overlay=False,background=False;
  int cross=False,circle=False,fill=False,inc;
  float *kptmap,fmax,radius=2.5;
  unsigned char *cmap,*keypoints,*ovlimage,*ovltmp;
  FILE *fd;

  /*==============================================================
   * Interpret parameters from Command-line Arguments 
   *==============================================================
   */
  
  if ((exists_argument(argc,argv,"-i"))){
    fgen=get_char_argument(argc,argv,"-i");
    
    if (exists_argument(argc,argv,"-th"))
      thresh=get_float_argument(argc,argv,"-th");
    if (exists_argument(argc,argv,"-size"))
      size=(int) get_float_argument(argc,argv,"-size");
        
    if (exists_argument(argc,argv,"-asc"))
      ascii=True;

    if (exists_argument(argc,argv,"-over")){
      overlay=True;
      type=get_char_argument(argc,argv,"-over");
    }

    if (exists_argument(argc,argv,"-back")) {
      fovl=get_char_argument(argc,argv,"-back");
      background=True;
    }

    if (exists_argument(argc,argv,"-rad"))
      radius=get_float_argument(argc,argv,"-rad");

  }
  else{
    printf(errmessage);
    exit(1);
  }
  
  /*================================================================
   * READ INPUT DATA (KEYMAP)
   *================================================================
   */
  
  sprintf(fout,"%s-kptmap.dat",fgen);
  if (fd=fopen(fout,"r")){
    fscanf(fd,"%s\n",arr);
    if (strcmp(arr,"BINARY")==0){
      fscanf(fd,"%d %d\n",&nx,&ny);
      kptmap=(float *) malloc(nx*ny*sizeof (float));
      fread((void *) kptmap,sizeof(float),nx*ny,fd);
      fprintf(stderr,"Loaded BINARY DATA (kptmap)\n",fout,nx,ny);
      fclose(fd);
    }
    fclose(fd);
  }
  else {
    sprintf(fout,"%s-kptmap.ras",fgen);
    if (read_image(&cmap,fout,&nx,&ny)){
      fprintf(stderr,"Loaded SUN RASTER DATA (kptmap)\n",fout,nx,ny);
      kptmap=(float *) malloc(nx*ny*sizeof (float));
      for (i=0;i<nx*ny;i++) kptmap[i]=(float) cmap[i];
      free(cmap);
    }
    else {
      fprintf(stderr,"Couldn't open file %s\n",fout);
      exit(1);
    }
  }

  if (overlay){
    if (background && read_image(&ovlimage,fovl,&nx2,&ny2))
      fprintf(stderr,"Loaded SUN RASTER IMAGE (overlay)\n",fovl,nx2,ny2);
    else {
      ovlimage=(unsigned char *) malloc(nx*ny*sizeof(char));
      memset((void *) ovlimage,0,nx*ny*sizeof(char));
      nx2=nx;ny2=ny;      
    }
    if ((nx!=nx2 || ny!=ny2)) overlay=False;
  }

  /*================================================================
   * THRESHOLDING AND LOCAL MAXIMA OPERATION
   *================================================================
   */

  s=(size-1)/2;
  fmax=-1e38;
  for (i=0;i<nx*ny;i++) if (kptmap[i] > fmax) fmax=kptmap[i];
  for (i=0;i<nx*ny;i++)
    if (kptmap[i] < fmax*thresh/100.0) kptmap[i]=0.0;
  
  keypoints=(unsigned char *) malloc(nx*ny*sizeof (char));
  memset((void *) keypoints,0,nx*ny*sizeof (char));
  
  for (y=s;y<ny-s;y++){
    for (x=s;x<nx-s;x++){
      if (kptmap[x+y*nx] > 0.0) {
	nms=True;
	for (dy=-s;dy<=s;dy++){
	  for (dx=-s;dx<=s;dx++)
	    if ((dx != 0 || dy != 0) &&
		kptmap[x+dx+(y+dy)*nx] >= kptmap[x+y*nx]) nms=False;
	}
	if (nms) keypoints[x+y*nx]=255;
      }
    }
  }

  /*================================================================
   * OUTPUT DATA
   *================================================================
   */

  sprintf(fout,"%s-kpts.ras",fgen);
  write_grey_image(keypoints,fout,nx,ny,sizeof(char),1);
  fprintf(stderr,"Saved sun-raster: %s\n",fout);
  
  if (overlay){
    
    ovltmp=(unsigned char *) malloc(nx*ny*sizeof (char));
    memset((void *) ovltmp,0,nx*ny*sizeof (char));
    d=IROUND(radius);
    
    if (strcmp(type,"cross")==0){ 
      for (y=d;y<(ny-d);y++){
	for (x=d;x<(nx-d);x++){
	  if (keypoints[x+y*nx]){
	    for (dx=-(d);dx<=(d);dx++)
	      ovltmp[x+dx+y*nx]=ovltmp[x+(y+dx)*nx]=keypoints[x+y*nx];
	  }
	}
      }
    }
    
    if (strcmp(type,"circle")==0){
      for (y=d;y<(ny-d);y++){
	for (x=d;x<(nx-d);x++){
	   if (keypoints[x+y*nx]){
	     for (inc=0;inc<=20;inc ++){
	       dx=(int) (radius*fcos(inc*0.05*M_PI));
	       dy=(int) (radius*fsin(inc*0.05*M_PI));
	       if ((fabs(dx) <= d)&&(fabs(dy) <= d))
		 ovltmp[x+dx+(y-dy)*nx]=ovltmp[x+dx+(y+dy)*nx]=keypoints[x+y*nx];
	     }
	   }
	}

      }
    }
    
    sprintf(fout,"%s-kpts-ovl.ras",fgen);
    write_overlay_image(ovlimage,ovltmp,fout,nx,ny,255,0,0);
    fprintf(stderr,"Saved OVERLAY sun-raster: %s\n",fout);
    free(ovltmp);free(ovlimage);
  }

  if (ascii){
    sprintf(fout,"%s-kpts.ascii",fgen);
    fd=fopen(fout,"w");
    for (y=0;y<ny;y++) {
      for (x=0;x<nx;x++) 
	if (keypoints[x+y*nx]) fprintf(fd,"%d %d %f\n",x,y,kptmap[x+y*nx]);
    }
    fclose(fd);
    fprintf(stderr,"Saved key-point ascii file: %s\n",fout);
  }

  /*================================================================
   * FREE ALLOCATED DATA
   *================================================================
   */

  free(keypoints);
  free(kptmap);

  exit(0);
}
