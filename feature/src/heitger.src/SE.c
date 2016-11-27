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

static char errmessage[] = "Usage:  SE -i image -o generic_filename [options]\n        ,where options are any of the following\n\n        -nori [int]       : number of filter orientations (def. 6)\n\n        -sigma [float]    : sigma of gaussian envelope (def. 2.0)\n        -sweep [float]    : sweep parameter of s-gabor filter (def. 0.50)\n        -orisel [int]     : orient. selectivity of gabor filter [2,4,6] (def. 2)\n\n        -version          : version (1 or 2) of SE operator (def. 2)\n        -supp [float]     : suppression factor for SE operator (def. 2.0)\n        -enh [float]      : enhancement factor for SE operator (def. 2.0)\n        -cosw [float]     : weight of even filter for SE operator (def. 1.0)\n\n        -keypoints        : compute key-points\n        -xgain [float]    : cross ridge gain for key-points (def. 2.5)\n        -xcorr [float]    : cross ridge correction for key-points (def. 0.18)\n        -compgain [float] : compensation gain for key-points (def. 1.0)\n        -evengain [float] : even gain for key-points (def. 2.0)\n\n        -quality          : compute general edge quality\n\n        -nobinary         : binary output is disabled\n        -raster           : all .dat output is also saved as SUN-raster\n\n        The program performs in sequence: filter the image with s-gabor filters,\n        compute oriented energy, apply the Suppression & Enhancement (SE) operator,\n        perform non-maximum suppression, and optionally compute key-point map,\n        and general edge quality.\n";
 

void main(int argc,char **argv)
{
  char *fin=NULL,fout[256],*fname=NULL;
  int i,nx,ny,tnx,tny;
  int version=2,raster=False,binary=True;
  int comp_quality=False,comp_keypoints=False;
  unsigned char *image;
  float sigma=2.0,sweep=0.5,orisel=2.0;   /* parameters energy filters   */
  int nori=6;                   /* number of filter orientations */
  float supp=2.0,enh=2.0,cosw=1.0;        /* parameters SE-operator */
  float xgain=2.50,xcorr=0.18,compgain=1.0,
    evengain=2.0;                 /* parameters key-point scheme */
  float *maxmap,*locori,*nms;
  unsigned char *type;
  float *kptmap,*quality,*mod_maxmap;
  FILE *fd;
  float *ftmp;

  /*==============================================================
   * Interpret parameters from Command-line Arguments 
   *==============================================================
   */

  
  if ((exists_argument(argc,argv,"-i"))&&(exists_argument(argc,argv,"-o"))){
    fin=get_char_argument(argc,argv,"-i");
    fname=get_char_argument(argc,argv,"-o");
    
    if (exists_argument(argc,argv,"-sigma"))
      sigma=get_float_argument(argc,argv,"-sig");
    if (exists_argument(argc,argv,"-sweep"))
      sweep=get_float_argument(argc,argv,"-sweep"); 
    if (exists_argument(argc,argv,"-orisel"))
      orisel=get_float_argument(argc,argv,"-orisel"); 
    
    if (exists_argument(argc,argv,"-nori"))
      nori=(int) get_float_argument(argc,argv,"-nori"); 
    
    if (exists_argument(argc,argv,"-supp"))
      supp=get_float_argument(argc,argv,"-supp"); 
    if (exists_argument(argc,argv,"-enh"))
      enh=get_float_argument(argc,argv,"-enh"); 
    if (exists_argument(argc,argv,"-cosw"))
      cosw=get_float_argument(argc,argv,"-cosw"); 

    if (exists_argument(argc,argv,"-xgain"))
      xgain=get_float_argument(argc,argv,"-xgain"); 
    if (exists_argument(argc,argv,"-xcorr"))
      xcorr=get_float_argument(argc,argv,"-xcorr"); 
    if (exists_argument(argc,argv,"-compgain"))
      compgain=get_float_argument(argc,argv,"-compgain"); 
    if (exists_argument(argc,argv,"-evengain"))
      evengain=get_float_argument(argc,argv,"-evengain");

    if (exists_argument(argc,argv,"-version"))
      version=(int) get_float_argument(argc,argv,"-version");
    
    if (exists_argument(argc,argv,"-nobin"))
      binary=False;
    if (exists_argument(argc,argv,"-ras"))
      raster=True;
    if (exists_argument(argc,argv,"-key"))
      comp_keypoints=True;
    if (exists_argument(argc,argv,"-qual"))
      comp_quality=True;
    
    

  }
  else{
    printf(errmessage);
    exit(1);
  }
  
  /*================================================================
   * READ SUN RASTER IMAGE
   *================================================================
   */

  if (read_image(&image,fin,&nx,&ny)){
    fprintf(stderr,"Loaded image: %s [%d,%d]\n",fin,nx,ny);
  }
  else {
    fprintf(stderr,"Couldn't open file %s\n",fin);
    exit(1);
  }

  /*================================================================
   * ALLOCATE MAPS
   *================================================================
   */
  
  maxmap=(float *) malloc(nx*ny*sizeof(float));
  locori=(float *) malloc(nx*ny*sizeof(float));
  if (comp_keypoints) kptmap=(float *) malloc(nx*ny*sizeof(float));
  else kptmap=NULL;
  if (comp_quality) quality=(float *) malloc(nx*ny*sizeof(float));
  else quality=NULL;
  mod_maxmap=(float *) malloc(nx*ny*sizeof(float));
  nms=(float *) malloc(nx*ny*sizeof (float));
  type=(unsigned char *) malloc(nx*ny*sizeof (char));
  
  /*================================================================
   * ECHO PARAMETERS
   *================================================================
   */
  
  fprintf(stderr,"\n======= PARAMETERS ======================\n");
  fprintf(stderr,"Version(%d): Image size [%d,%d]\n",version,nx,ny);
  fprintf(stderr,"Number of filter orientations: %d\n",nori);
  fprintf(stderr,"Egabor  : sigma=%.2f,sweep=%.2f,orisel=%.0f\n",
	  sigma,sweep,orisel);
  fprintf(stderr,"S&E oper: supp=%.1f,enh=%.1f,cosw=%.2f\n",
	  supp,enh,cosw);
  fprintf(stderr,"Keypoint: xgain=%.2f,xcorr=%.2f\n",xgain,xcorr);
  fprintf(stderr,"Keypoint: comp-gain=%.1f,even-gain=%.1f\n",
	  compgain,evengain);
  fprintf(stderr,"=========================================\n");
  
  /*================================================================
   * COMPUTE FEATURE EXTRACTION
   *================================================================
   */

  Low_level_Processing(image,nx,ny,nori,
		       sigma,sweep,orisel,
		       supp,enh,cosw,
		       xgain,xcorr,compgain,evengain,
		       maxmap,mod_maxmap,locori,nms,
		       type,kptmap,quality,(version==2)?2:1);
     
  /*================================================================
   * OUTPUT DATA 
   *================================================================
   */

  if (raster) {
    sprintf(fout,"%s-nms.ras",fname);
    write_grey_image(nms,fout,nx,ny,sizeof(float),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
    sprintf(fout,"%s-maxmap.ras",fname);
    write_grey_image(maxmap,fout,nx,ny,sizeof(float),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
    sprintf(fout,"%s-mod_maxmap.ras",fname);
    write_grey_image(mod_maxmap,fout,nx,ny,sizeof(float),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
    sprintf(fout,"%s-locori.ras",fname);
    write_grey_image(locori,fout,nx,ny,sizeof(float),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
    if (comp_keypoints){
      sprintf(fout,"%s-kptmap.ras",fname);
      write_grey_image(kptmap,fout,nx,ny,sizeof(float),1);
      fprintf(stderr,"Saved sun-raster: %s\n",fout);
    }
    if (comp_quality){
      sprintf(fout,"%s-quality.ras",fname);
      write_grey_image(quality,fout,nx,ny,sizeof(float),1);
      fprintf(stderr,"Saved sun-raster: %s\n",fout);
    }
    sprintf(fout,"%s-type.ras",fname);
    write_grey_image(type,fout,nx,ny,sizeof(char),1);
    fprintf(stderr,"Saved sun-raster: %s\n",fout);
  }

  if (binary) {
    sprintf(fout,"%s-nms.dat",fname);
    fd=fopen(fout,"w");
    fprintf(fd,"BINARY\n%d %d\n",nx,ny);
    fwrite((void *) nms,sizeof(float),nx*ny,fd);
    fclose(fd);
    fprintf(stderr,"Saved binary data: %s\n",fout);

    sprintf(fout,"%s-maxmap.dat",fname);
    fd=fopen(fout,"w");
    fprintf(fd,"BINARY\n%d %d\n",nx,ny);
    fwrite((void *) maxmap,sizeof(float),nx*ny,fd);
    fclose(fd);
    fprintf(stderr,"Saved binary data: %s\n",fout);

    sprintf(fout,"%s-mod_maxmap.dat",fname);
    fd=fopen(fout,"w");
    fprintf(fd,"BINARY\n%d %d\n",nx,ny);
    fwrite((void *) mod_maxmap,sizeof(float),nx*ny,fd);
    fclose(fd);
    fprintf(stderr,"Saved binary data: %s\n",fout);

    sprintf(fout,"%s-locori.dat",fname);
    fd=fopen(fout,"w");
    fprintf(fd,"BINARY\n%d %d\n",nx,ny);
    fwrite((void *) locori,sizeof(float),nx*ny,fd);
    fclose(fd);
    fprintf(stderr,"Saved binary data: %s\n",fout);

    if (comp_keypoints){
      sprintf(fout,"%s-kptmap.dat",fname);
      fd=fopen(fout,"w");
      fprintf(fd,"BINARY\n%d %d\n",nx,ny);
      fwrite((void *) kptmap,sizeof(float),nx*ny,fd);
      fclose(fd);
      fprintf(stderr,"Saved binary data: %s\n",fout);
    }
    
    if (comp_quality){
      sprintf(fout,"%s-quality.dat",fname);
      fd=fopen(fout,"w");
      fprintf(fd,"BINARY\n%d %d\n",nx,ny);
      fwrite((void *) quality,sizeof(float),nx*ny,fd);
      fclose(fd);
      fprintf(stderr,"Saved binary data: %s\n",fout);
    }

    if (!raster) {
      sprintf(fout,"%s-type.ras",fname);
      write_grey_image(type,fout,nx,ny,sizeof(char),1);
      fprintf(stderr,"Saved sun-raster: %s\n",fout);
    }
  }

  /*================================================================
   * FREE ALLOCATED DATA
   *================================================================
   */

  free(maxmap);
  free(locori);
  if (kptmap) free(kptmap);
  if (quality) free(quality);
  free(mod_maxmap);
  free(nms);
  free(type);
  
  exit(0);
}
