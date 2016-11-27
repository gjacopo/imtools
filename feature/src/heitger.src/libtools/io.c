/* 
 * $Log: io.c,v $
 * Revision 1.3  1995/11/09  12:56:58  ohenri
 * *** empty log message ***
 *
 * Revision 1.2  1995/10/02  12:25:57  ohenri
 * *** empty log message ***
 *
 * Revision 1.1  1995/06/21  12:58:19  ohenri
 * Initial revision
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
#include <malloc.h>
#include "tools.h"


/*****************************************************************************
 * int read_image(unsigned char **image,char fname[],int *nx, int *ny)
 *----------------------------------------------------------------------------
 * read_image() function reads an image saved as a SUN rasterfile.
 * Call this routine with an unsigned character pointer and the program 
 * returns the allocated image, with its size (nx,ny).
 *****************************************************************************
 */

int read_image(unsigned char **image,char *fname,int *nx,int *ny)
{
  FILE  *file_id;
  int header[8],i,i_size,ret=0,dx,dy,x,y;
  unsigned char *colormap=NULL,*R,*G,*B,tmp;
 
  file_id=fopen(fname,"r");
  if (file_id == NULL){
    printf("Error: Can't open file: %s\n",fname);
    return 0;;
  }
  else{
    fread(header,sizeof(int),8,file_id);
    dx=*nx=header[1];
    dy=*ny=header[2];
    i_size=(*nx)*(*ny);
    
    if (header[5] == RT_STANDARD){
      
      *image=(unsigned char *) malloc((unsigned) i_size*sizeof(char));
      switch (header[3]){
      case 8:
	if ((header[6] == RMT_EQUAL_RGB)&&(header[7])){
	  colormap=(unsigned char *) malloc(header[7]*sizeof(char));
	  fread(colormap,1,header[7],file_id);
	  /* fprintf(stderr,"Mapping Colormap (%d bytes).\n",header[7]); */
	}
	
	for (i=0;i<*ny;i++){
	  fread((*image+i*(*nx)),1,(*nx),file_id);
	  if ((*nx)%2) fread(&tmp,1,1,file_id);
	}
	
	if (colormap) 
	  for (i=0;i<i_size;i++) (*image)[i]=colormap[(*image)[i]];
	
	ret=1;break;
      case 24:
	if (header[6] == RMT_NONE){
	  
	  R=(unsigned char *) malloc(dx*dy*sizeof(char));
	  G=(unsigned char *) malloc(dx*dy*sizeof(char));
	  B=(unsigned char *) malloc(dx*dy*sizeof(char));
      
	  for (y=0;y<dy;y++){
	    for (x=0;x<dx;x++){
	      fread(&(B[x+y*dx]),1,1,file_id);
	      fread(&(G[x+y*dx]),1,1,file_id);
	      fread(&(R[x+y*dx]),1,1,file_id);
	    }
	    
	    if (dx%2) {
	      fread(&tmp,1,1,file_id);
	    }
	  }

	  /* Value channel is used as grey-value */
	  for (i=0;i<dx*dy;i++) 
	    (*image)[i]=(unsigned char) MAX((int) MAX(R[i],G[i]),(int) B[i]);
	  free(R);free(G);free(B);
	}
	ret=1;
	break;
      default: fprintf(stderr,
      "Error reading Sun raster file: unsupported depth=%d bit.\n",header[3]);
	return 0;
      }
      
      fclose(file_id);
      
      return ret;
    }
    else{
      fprintf(stderr,"Error: not a standard SUN raster file.\n");
      return 0;
    }
  }
}

int read_RGB_image(unsigned char **R,unsigned char **G,unsigned char **B,
		   char *fname,int *nx,int *ny)
{
  FILE  *file_id;
  int header[8],i,i_size,x,y,dx,dy;
  unsigned char tmp;
  
  file_id=fopen(fname,"r");
  if (file_id == NULL){
    printf("Error: Can't open file: %s\n",fname);
    return 0;
  }
  else{
    fread(header,sizeof(int),8,file_id);
    dx=*nx=header[1];
    dy=*ny=header[2];
    
    if (header[5] == RT_STANDARD && header[3] == 24){
      
      *R=(unsigned char *) malloc(dx*dy*sizeof(char));
      *G=(unsigned char *) malloc(dx*dy*sizeof(char));
      *B=(unsigned char *) malloc(dx*dy*sizeof(char));
      
      for (y=0;y<dy;y++){
	for (x=0;x<dx;x++){
	  fread((*B+x+y*dx),1,1,file_id);
	  fread((*G+x+y*dx),1,1,file_id);
	  fread((*R+x+y*dx),1,1,file_id);
	}
	
	if (dx%2) {
	  fread(&tmp,1,1,file_id);
	}
      }
      fclose(file_id);
      
      return 1;
    }
    else{
      fprintf(stderr,"Error: not a 24 bits SUN raster file.\n");
      return 0;
    }
  }
}

 
int write_grey_image(void *image,char *filename,int nx,int ny,
		      int type,int colormap)
{
  FILE *file_id;
  float *fimg,mmax=-1e38,scale=0.0,mmin=1e38;
  unsigned char *cimg,tmp=0,cmap[256];
  int i,width;
  static int Header[] = { RAS_MAGIC, 0, 0, 8, 0, RT_STANDARD, RMT_NONE, 0};
  
  Header[1]=nx;
  Header[2]=ny;
  Header[4]=nx*ny;
  if (colormap){
    Header[6]=RMT_EQUAL_RGB;
    Header[7]=3*256;
    for (i=0;i<256;i++) cmap[i]=i;
  }
  
  if (type == FLOAT){ /* Scale image to 8 bit representation */
    cimg=(unsigned char *) malloc((unsigned) Header[4]*sizeof(char));
    fimg=(float *) image;
    for (i=0;i<Header[4];i++){
      if (fimg[i] > mmax) mmax=fimg[i];
      if (fimg[i] < mmin) mmin=fimg[i];
    }
    if ((mmax-mmin) !=  0.0) scale=255.0/(mmax-mmin);
    for (i=0;i<Header[4];i++) cimg[i]=(unsigned char) ((fimg[i]-mmin)*scale);
  }
  else if (type == CHAR) cimg=(unsigned char *) image;
  else {fprintf(stderr,"write_grey_image(): type=%d error\n",type);return 0;}
  
  if ((file_id=fopen(filename,"w")) != NULL){
    fwrite( Header, sizeof(int), 8, file_id);
    if (colormap) {for (i=0;i<3;i++) fwrite(cmap,256,1,file_id);}
    for (i=0;i<ny;i++){
      fwrite((cimg+i*nx),1,nx,file_id);
      if (nx%2) fwrite(&tmp,1,1,file_id);
    }
    fclose(file_id);
  }
  else {fprintf(stderr,"write_grey_image(): can't open file\n");return 0;}
  
  if (type == FLOAT) free(cimg);
  return 1;
}

int write_RGB_image(unsigned char *R,unsigned char *G,unsigned char *B,
		    char *filename,int nx,int ny)
{
  FILE *file_id;
  int x,y;
  static int Header[] = { RAS_MAGIC, 0, 0, 24, 0, RT_STANDARD, RMT_NONE, 0};
  unsigned char tmp=0;

  Header[1]=nx;
  Header[2]=ny;
  Header[4]=3*nx*ny;
  
  if ((file_id=fopen(filename,"w")) != NULL){
    fwrite(Header, sizeof(int), 8, file_id);
    for (y=0;y<ny;y++){
      for (x=0;x<nx;x++){
	fwrite(&(B[x+y*nx]),1,1,file_id);
	fwrite(&(G[x+y*nx]),1,1,file_id);
	fwrite(&(R[x+y*nx]),1,1,file_id);
      }
      
      if (nx%2) {
	fwrite(&tmp,1,1,file_id);
      }
    }
    fclose(file_id);
  }
  else {fprintf(stderr,"write_RGB_image(): can't open file\n");return 0;}
  
  return 1;
}


void write_color_image(unsigned char *image,char *fname,int nx,int ny)
{
  FILE *file_id;
  static int Header[]={RAS_MAGIC, 0, 0, 8, 0, RT_STANDARD, RMT_EQUAL_RGB, 0};
  int width,i;
  char tmp=0;
  static unsigned char 
    red[NCMAP]={0, 80,140,255,255,  0,  0,160,240,255,160},
    green[NCMAP]={0,170,140,160,255,255,  0,  0,  0,255,120},
    blue[NCMAP]={0,255,140,  0,  0,  0,255,200,  0,255,  0};
  /* [black,lightblue,grey,orange,yellow,green,blue,violett,red,white,brown] */
  
  file_id=fopen(fname,"w");
  
  Header[1]=nx;
  Header[2]=ny;
  Header[4]=nx*ny;
  Header[7]=3*NCMAP;

  fwrite( Header, sizeof(int), 8, file_id);
  fwrite(red  ,NCMAP,1,file_id);
  fwrite(green,NCMAP,1,file_id);
  fwrite(blue ,NCMAP,1,file_id);

  for (i=0;i<ny;i++){
    fwrite((image+i*nx),1,nx,file_id);
    if (nx%2) fwrite(&tmp,1,1,file_id);
  }
  fclose(file_id);
}

void write_overlay_image(unsigned char *image,unsigned char *overlay,
			 char *fname,int nx,int ny,char R,char G,char B)
{
  int i,i_size=nx*ny,width;
  FILE *file_id;
  static int Header[8]={RAS_MAGIC, 0, 0, 8, 0, RT_STANDARD, RMT_EQUAL_RGB, 0};
  unsigned char *tmp,l,ctmp=0;
  unsigned char red[256],green[256],blue[256];

  tmp=(unsigned char *) malloc(i_size*sizeof(unsigned char));
  for (i=0;i<i_size;i++){ 
    if (image[i] == 255) tmp[i]=254;
    else tmp[i]=(unsigned char) image[i];
    if (overlay[i] > 0) tmp[i]=255; 
  }    

  /* make greyscale colormap of image */

  for (l=0;l<255;l++){
    red[l]=green[l]=blue[l]=l;
  }
  red[255]=R;green[255]=G;blue[255]=B;

  file_id=fopen(fname,"w");
  
  Header[1]=nx;
  Header[2]=ny;
  Header[4]=i_size;
  Header[7]=3*256;

  fwrite( Header, sizeof(int), 8, file_id);
  fwrite(red  ,256,1,file_id);
  fwrite(green,256,1,file_id);
  fwrite(blue ,256,1,file_id);
  
  for (i=0;i<ny;i++){
    fwrite((tmp+i*nx),1,nx,file_id);
    if (nx%2) fwrite(&ctmp,1,1,file_id);
  }
  free(tmp);
  fclose(file_id);
}

char *get_generic(char *string)
{
  int i,len;
  char *ret;

  len=strlen(string);
  ret=(char *) malloc(len*sizeof(char));
  memcpy((void *) ret,(void *) string,len*sizeof(char));

  if (ret){
    for(i=len;i>=0;i--)
      if (ret[i] == '.'){
	ret[i]='\0';
	return ret;
      }
    return ret;
  }

  return (char *) NULL;
}


int exists_argument(int argc,char **argv,char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strstr(argv[i],keystr)) return 1;
  return 0;
}

float get_float_argument(int argc,char **argv,char *keystr)
{
  int i;
  float result=0.0;

  for (i=1;i<argc;i++){
    if (strstr(argv[i],keystr)){
      result=atof(argv[i+1]);
      break;
    }
  }
  if (i == argc){ 
    fprintf(stderr,"Error: get_float_argument(%s)\n",keystr);
    return (float) -999999999;
  }
  else return result;
}

char *get_char_argument(int argc,char **argv,char *keystr)
{
  int i;
  char *result;

  result=(char *) malloc(256);
  result[0]=0;

  for (i=1;i<argc;i++){
    if (strstr(argv[i],keystr)){
      sprintf(result,"%s",argv[i+1]);
      break;
    }
  }
  if (i == argc) fprintf(stderr,"Error: get_char_argument(%s)\n",keystr);

  return result;
}

