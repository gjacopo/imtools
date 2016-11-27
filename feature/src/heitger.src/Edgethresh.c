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

static char errmessage[] = "Usage:  Edgethresh -i generic_filename [options]\n        ,where options are any of the following\n\n        -sign [float]        : significance level (in %%) (def. 99.0%%)\n        -fixed [float]       : threshold (in %%) relative to maximum (def. 10.0%%)\n                               (this option overrides -sign)\n\n        -ascii               : output is an ascii file with [x,y,magnitude]\n                               if type-input exists [x,y,magnitude,type]\n\n        -overlay [filename]  : output is an (color) overlaid raster file\n\n        The program estimates the image noise, based on the method of\n        Voorhees and Poggio, IU Workshop 1987. Given the significance\n        level (def. 99%%), an edge detection threshold level is computed.\n        The generic filename must be the same as the one used after -o\n        for the SE program!\n        \n";

void main(int argc,char **argv)
{
  char fout[256],*fovl,arr[256],*fgen,*type;
  float signif=99.0;
  int nx,ny,nx2,ny2,i,x,y;
  int ascii=False,overlay=False;
  int fixed=False,typecomp=True;
  float *fedgemap,*fmaxmap,*fnms,*ftype;
  float mode,dmax,dmin,thresh,fthresh;
  unsigned char *edgemap,*cnms,*cmaxmap,*ctype=NULL;
  unsigned char *ovlimage,*ovltmp;
  FILE *fd;

  /*==============================================================
   * Interpret parameters from Command-line Arguments 
   *==============================================================
   */
  
  if (exists_argument(argc,argv,"-i")){
      fgen=get_char_argument(argc,argv,"-i");

    if (exists_argument(argc,argv,"-sign"))
      signif=get_float_argument(argc,argv,"-sign");
    if (exists_argument(argc,argv,"-fixed")){
      fthresh=get_float_argument(argc,argv,"-fixed");
      fixed=True;
    }
    
    if (exists_argument(argc,argv,"-asc"))
      ascii=True;
    if (exists_argument(argc,argv,"-overlay")){
      overlay=True;
      fovl=get_char_argument(argc,argv,"-overlay");
    }
  }
  else{
    printf(errmessage);
    exit(1);
  }
  
  /*================================================================
   * READ INPUT DATA (NMS and MAXMAP)
   *================================================================
   */
  
  /* read NMS map */
  
  sprintf(fout,"%s-nms.dat",fgen);
  if (fd=fopen(fout,"r")){
    fscanf(fd,"%s\n",arr);
    if (strcmp(arr,"BINARY")==0){
      fscanf(fd,"%d %d\n",&nx,&ny);
      fnms=(float *) malloc(nx*ny*sizeof (float));
      fread((void *) fnms,sizeof(float),nx*ny,fd);
      fprintf(stderr,"Loaded BINARY DATA (nms)\n",fout,nx,ny);
      fclose(fd);
    }
    else {
      fprintf(stderr,"Wrong format of BINARY file %s\n",fout);
      exit(1);
    }
  }
  else {
    fprintf(stderr,"Couldn't open BINARY file %s\n",fout);
    exit(1);
  }
  
  /* read MAXMAP map */
  
  sprintf(fout,"%s-maxmap.dat",fgen);
  if (fd=fopen(fout,"r")){
    fscanf(fd,"%s\n",arr);
    if (strcmp(arr,"BINARY")==0){
      fscanf(fd,"%d %d\n",&nx,&ny);
      fmaxmap=(float *) malloc(nx*ny*sizeof (float));
      fread((void *) fmaxmap,sizeof(float),nx*ny,fd);
      fprintf(stderr,"Loaded BINARY DATA (maxmap)\n",fout,nx,ny);
      fclose(fd);
    }
    else {
      fprintf(stderr,"Wrong format of BINARY file %s\n",fout);
      exit(1);
    }
  }
  else {
    fprintf(stderr,"Couldn't open BINARY file %s\n",fout);
    exit(1);
  }

  /* read optional TYPE map */
  
  sprintf(fout,"%s-type.ras",fgen);
  if (read_image(&ctype,fout,&nx,&ny))
    fprintf(stderr,"Loaded SUN RASTER DATA (type)\n",fout,nx,ny);
  else typecomp=False;
  
  /* read optional OVERLAY map */
  
  if (overlay){
    if (read_image(&ovlimage,fovl,&nx2,&ny2))
      fprintf(stderr,"Loaded SUN RASTER IMAGE (overlay)\n",fovl,nx2,ny2);
    else {
      ovlimage=(unsigned char *) malloc(nx*ny*sizeof(char));
      memset((void *) ovlimage,0,nx*ny*sizeof(char));
      nx2=nx;ny2=ny;      
    }
    if ((nx!=nx2 || ny!=ny2)) overlay=False;
  }

  /*================================================================
   * NOISE ESTIMATION AND THRESHOLDING
   *================================================================
   */
  
  edgemap=(unsigned char *) malloc(nx*ny*sizeof(char));
  memset((void *) edgemap,0,nx*ny*sizeof(char));
  
  NoiseEstimation(fmaxmap,nx,ny,signif,False,&mode,&dmax,&dmin,&thresh);
    
  if (!fixed)
    fprintf(stderr,"Min/Max Value: [%.3f,%.3f]\nMode=%f (%.1f %%)\nSignificance factor: %.1f\nAutomatic threshold: %.1f %%\n",dmin,dmax,mode,mode/dmax*100.0,fsqrt(-2.0*log((100.0-signif)/100.0)),thresh/dmax*100.0);
  else fprintf(stderr,"Min/Max Value: [%.3f,%.3f]\nFixed threshold: %.1f %%\n",dmin,dmax,fthresh);
    
  for (i=0;i<nx*ny;i++) { 
    if (!fixed && fnms[i] >= thresh) edgemap[i]=255;
    else if (fixed && fnms[i] >= fthresh*dmax/100.0) edgemap[i]=255;
  }

  /*================================================================
   * OUTPUT DATA
   *================================================================
   */

  if (typecomp){
    for (i=0;i<nx*ny;i++) if (edgemap[i]) edgemap[i]=ctype[i];
    sprintf(fout,"%s-edgemap.ras",fgen);
    write_grey_image(edgemap,fout,nx,ny,sizeof(char),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
  }
  else {
    sprintf(fout,"%s-edgemap.ras",fgen);
    write_grey_image(edgemap,fout,nx,ny,sizeof(char),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
  }
   
  if (overlay){
    sprintf(fout,"%s-edgemap-ovl.ras",fgen);
    write_overlay_image(ovlimage,edgemap,fout,nx,ny,255,0,0);
    fprintf(stderr,"Saved OVERLAY sun-raster: %s\n",fout);
    free(ovlimage);
  }
  
  if (ascii){
    sprintf(fout,"%s-edgemap.ascii",fgen);
    fd=fopen(fout,"w");
    for (y=0;y<ny;y++) {
      for (x=0;x<nx;x++) {
	if (edgemap[x+y*nx] && typecomp) 
	  fprintf(fd,"%d %d %f %d\n",x,y,fmaxmap[x+y*nx],
		  (ctype[x+y*nx]==EDGETYPE)?200:100);
	else if (edgemap[x+y*nx]) 
	  fprintf(fd,"%d %d %f 0\n",x,y,fmaxmap[x+y*nx]);
      }
    }
    fclose(fd);
    fprintf(stderr,"Saved edge pixels ascii file: %s\n",fout);
  }
  
  /*================================================================
   * FREE ALLOCATED DATA
   *================================================================
   */

  free(fnms);
  free(fmaxmap);
  free(edgemap);
  if (typecomp) free(ctype);
  
  exit(0);
}
