/* has the code to implement the optimal zero crossing edge detector */

/************************************************************************/
#include "header.h"

static float *temp_array, *temp_array1, *temp_array2;
float canny_op[500], canny_sm[500];
int W;
int canny_mask (float sigma, int mrc);

int zero_cross_filter() 
{
  float    *derivx_ptr, *derivy_ptr, *temp_ptr, *dest, *source1, *source2;
  float     temp1, temp2, w_1, w_2, max_val;
  UCHAR    *dest1;
  char file2[100] ;
  register int       i, j;

  /*fprintf(stderr, "\n OPERATOR : %d", OPERATOR);*/

  if (OPERATOR == 2) {
    canny_mask (Gamma, MRC-1);
  }
  if ((temp_array = (float *)calloc((max(expnd_width, expnd_height)),
				     sizeof(float))) == NULL) {
    fprintf(stderr,"\n Problems allocating memory for a temporary vector");
      exit(1) ;  
    }

  if ((temp_array1 = (float *)calloc((max(expnd_width, expnd_height)),
				     sizeof(float))) == NULL) {
    fprintf(stderr,"\n Problems allocating memory for a temporary vector");
    exit(1) ;  
  }
  if ((temp_array2 = (float *)calloc((max(expnd_width, expnd_height)),
				     sizeof(float))) == NULL) {
    fprintf(stderr,"\n Problems allocating memory for a temporary vector");
    exit(1) ;  
  }

  if ((derivy_ptr = (float *)calloc((exp_img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for filtered image") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  /* Do filtering to evaluate derivative in the y (row) direction */
  
  if(derivy_convl(image_ptr, derivy_ptr) == NOOK)
    {
      fprintf(stderr, "Failure in row recursive filtering.\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  

  if ((derivx_ptr = (float *)calloc((exp_img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for filtered image") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  /* Do filtering to evaluate derivative in the x (column) direction */
  
  
  if(derivx_convl(image_ptr, derivx_ptr) == NOOK)
    {
      fprintf(stderr, "Failure in column recursive filtering.\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  
  free(temp_array); free(temp_array2);  free(temp_array1);
  
    
  if ((temp_ptr = (float *)calloc((exp_img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for a temporary array \n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  
  /*save_convl (derivx_ptr, "tempx.out");
    save_convl (derivy_ptr, "tempy.out");*/

  dest = temp_ptr;
  source1 = derivx_ptr;
  source2 = derivy_ptr; max_val = 0.0;
  for( i = 0 ; i < exp_img_size ; i++)
    {
      *dest++ = ( *source1 ) * ( *source1 ) + ( *source2 ) * ( *source2 );
      source1++;
      source2++;
    }

  dest = temp_ptr + expnd_width + 1;
  source1 = derivx_ptr + expnd_width + 1;
  source2 = derivy_ptr + expnd_width + 1;
  norm_const = 2* (3 * alpha + 5)/(beta*beta*beta);
  for ( i = expnd_width + 1 ; i < ( exp_img_size - expnd_width ) ; i++ )
    {
      temp1=  ( *source2 * tau_x  * ( *( dest + expnd_width) 
				      - *( dest - expnd_width))  
		+ *source1 * tau_y * ( *( dest + 1) - *( dest - 1))) 
	/ (2 * tau_x * tau_x  * *dest );
      
      /* correct for rounding off errors in regions of constant 
	 brightness */

      if (fabs(*dest) < 0.001) temp1 = 0.0; 

      temp1 = temp1/norm_const;
      if ((fabs(*source2) < 0.001)&& (fabs(*source1) < 0.001))
	temp2 = 100.1;
      else 
	temp2 = atan2(-*source2, -*source1);
      
      *source1++ = temp1;
      *source2++ = temp2;
      dest++;
    }
  /*save_convl (derivx_ptr, "temp.out");
    save_convl (derivy_ptr, "tempd.out");*/

  strcpy(file2, convlfilename) ;
  strcat(file2, ".num") ;

  /* *derivx_ptr now has zero crossing convolved image and 
     *derivy_ptr has the tangent direction in radians */

  free(temp_ptr);
  display_main ("Detecting zero crossings and authenticating them.");
  if (!(edge_ptr == NULL)) {free (edge_ptr);}
  if ((edge_ptr = (UCHAR *)calloc((UINT)(img_size), sizeof(UCHAR))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for edge array ") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  dest1=edge_ptr;
  for (i=0; i<img_size; i++)
    { *(dest1 +i) = 0;}
  
  if ((slope_ptr = (float *)calloc((img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for filtered image") ;
      exit(1) ;  
    }
  dest = slope_ptr; max_val = 0.0;
/**** The directions are set as follows 
         7(NW) | 0(N) | 1(NE)
	 -------------------
	 6(W)  |  *   | 2(E)
	 ------------------
	 5(SW) | 4(S) | 3(SE)
 The above denote the next pixel direction, so that the brighter
region is on the right ********/
  
  for( i=0; i<img_height; i++)
    { for( j=0; j<img_width; j++, dest1++, dest++)
      { source1 = derivx_ptr + (i + gksize) * expnd_width + j + gksize;
      source2 = derivy_ptr + (i + gksize) * expnd_width + j + gksize;
	  /* *derivx_ptr now has zero crossing convolved image and 
             *derivy_ptr has the tangent direction in radians */

	  if ((*source2 < 100.0) && (*source1 >= 0) && 
	      ((*(source1 - 1) < 0)||(*(source1 +1) < 0)
	       ||(*(source1 - expnd_width)< 0)
	       ||(*(source1 + expnd_width) < 0)))
	    { /* *dest1 = 255;*/
	      /*  fprintf (stderr,"%f ", *source2);*/
	       if (*source2 >= 0)
		{ if (*source2 < (pi/2))
		    {if (*source2 >= (pi/4))
		       
		       { w_1 =  cot(*source2)* (tau_y/tau_x);
			 w_2 =  (1 - w_1);
			 temp1 = (*(source1 + 1 + expnd_width) * w_1
				  + *(source1 + expnd_width) * w_2);
			 temp2 = (*(source1 - 1 - expnd_width) * w_1
				  + *(source1 - expnd_width) * w_2 ) ;
			if (((fabs(*source1) > 0.0001) || (temp1 > 0.0)) &&
			    (temp2 < 0))
			  {
			    if (fabs (*source2) < 3*pi/8)
			      /** between pi/4 and 3pi/8 **/
			      {*dest1 = 241;}
			    else
			      /** between 3pi/8 and pi/2 **/
			      {*dest1 = 242;}
			    *dest = (*source1 - temp2)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			    
			  }
		      }
		     else
		       { w_1 =  tan(*source2) * (tau_x/ tau_y);
			 w_2 = 1 - w_1;
			 temp1 = (*(source1 + 1 + expnd_width) * w_1 
				  + *(source1 +1) * w_2); 
			 temp2 = (*(source1 - 1 - expnd_width) * w_1
				  + *(source1 - 1) * w_2 );
			if (((fabs(*source1) > 0.0001) || (temp1 > 0.0)) &&
			    (temp2 < 0))
			  {
			    if (fabs (*source2) < pi/8)
			      /** between 0 and pi/8 **/
			      {*dest1 = 240;}
			    else
			      /** between pi/8 and pi/4 **/
			      {*dest1 = 241;}
			    *dest = (*source1 - temp2)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			  }
		      }
		   }
		else
		    {if (*source2 < (3 * pi/4))
		       { w_1 =  cot(*source2) *(- tau_y)/tau_x;
			 w_2 = 1- w_1;
			 temp1 = (*(source1 - 1 + expnd_width) * w_1
				 + *(source1 + expnd_width) * w_2); 
			 temp2 = (*(source1 + 1 - expnd_width) * w_1
				 + *(source1 - expnd_width) * w_2);
			if (((fabs(*source1) > 0.0001) || (temp1 > 0.0)) &&
			    (temp2 < 0))
			  {
			    if (fabs (*source2) < 5*pi/8)
			      /** between pi/2 and 5pi/8 **/
			      {*dest1 = 242;}
			    else
			      /** between 5pi/8 and 3pi/4 **/
			      {*dest1 = 243;}
			    *dest = (*source1 - temp2)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			  }
		      }
		     else
		       { w_1 = tan(*source2) * (- tau_x) / tau_y;
			 w_2 = 1 - w_1;
			 temp1 = (*(source1 - 1 + expnd_width) * w_1
				  + *(source1 - 1) * w_2);
			 temp2 = (*(source1 + 1 - expnd_width) * w_1
				  + *(source1 + 1) * w_2);
			if (((fabs(*source1) > 0.0001) || (temp1 > 0.0)) &&
			    (temp2 < 0))
			  {
			    if (fabs (*source2) < 7 * pi/8)
			      /** between 3pi/4 and 7pi/8 **/
			      {*dest1 = 243;}
			    else
			      /** between 7pi/8 and pi **/
			      {*dest1 = 244;}
			    *dest = (*source1 - temp2)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			  
			  }
		      }
		   }
		}
	      else 
		{ if (*source2 > -(pi/2))
		    {if (*source2 > - (pi/4))
		       { w_1 = tan(*source2) * (- tau_x) / tau_y;
			 w_2 = 1 - w_1;
			 temp1 = (*(source1 - 1 + expnd_width) * w_1
				  +  *(source1 - 1) * w_2);
			 temp2 = (*(source1 + 1 - expnd_width) * w_1
				  + *(source1 + 1) * w_2);
			if (((fabs(*source1) > 0.0001) || (temp2 > 0.0)) &&
			    (temp1 < 0))
			  {
			    if (fabs (*source2) < pi/8)
			      /** between 0 and -pi/8 **/
			      {*dest1 = 240;}
			    else
			      /** between -pi/8 and -pi/4 **/
			      {*dest1 = 247;}
			    *dest = (*source1 - temp1)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			  }
		      }
		     
		    else
		      { w_1 =  cot(*source2) *(- tau_y)/tau_x;
			w_2 = 1- w_1;
			temp1 = (*(source1 - 1 + expnd_width) * w_1
				 + *(source1 + expnd_width) * w_2);
			temp2 = (*(source1 + 1 - expnd_width) * w_1
				 + *(source1 - expnd_width) * w_2);
			if (((fabs(*source1) > 0.0001) ||  (temp2 > 0.0)) &&
			    (temp1 < 0))
			  {
			    if (fabs (*source2) < 3*pi/8)
			      /** between -pi/4 and -3pi/8 **/
			      {*dest1 = 247;}
			    else
			      {*dest1 = 246;}
			    /** between -3pi/8 and -pi/2 **/
			    *dest = (*source1 - temp1)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			  }
		      }
		   }   
		else
		  {if (*source2 > - (3 * pi/4))
		     { w_1 =  cot(*source2)* (tau_y/tau_x);
		       w_2 =  (1 - w_1);
		       temp1 = (*(source1 + 1 + expnd_width) * w_1
				+ *(source1 + expnd_width) * w_2);
		       temp2 = (*(source1 - 1 - expnd_width) * w_1
				+ *(source1 - expnd_width) * w_2) ;
		       if (((fabs(*source1) > 0.0001) || (temp2 > 0.0)) &&
			   (temp1 < 0))
			 {
			   if (fabs (*source2) < 5*pi/8)
			     /** between -pi/2 and -5pi/8 **/
			     {*dest1 = 246;}
			   else
			     /** between -5pi/8 and -3pi/4 **/
			     {*dest1 = 245;}
			   *dest = (*source1 - temp1)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			 }
		     }
		   
		  else
		    { w_1 =  tan(*source2) * (tau_x/ tau_y);
		      w_2 = 1 - w_1;
		      temp1 = (*(source1 + 1 + expnd_width) * w_1 
			       + *(source1 +1) * w_2);
		      temp2 = (*(source1 - 1 - expnd_width) * w_1
			       + *(source1 - 1) * w_2);
		      if (((fabs(*source1) > 0.0001) || (temp2 > 0.0)) &&
			  (temp1 < 0))
			{
			  if (fabs (*source2) < 7 * pi/8)
			    /** between -3pi/4 and -7pi/8 **/
			    {*dest1 = 245;}
			  else
			    /** between -7/4 and -pi **/
			    {*dest1 = 244;}
			  *dest = (*source1 - temp1)
			      /sqrt(w_1*w_1*tau_x*tau_x + tau_y * tau_y);
			}
		    }
		 }
		}
	       
	       if (max_val < fabs(*dest)) max_val = fabs(*dest);
	       /* fprintf (stderr,"%d \n", *dest1);
		  if ((fabs(w_1) > 1.01) || (fabs(w_2) > 1.01) ||
		  (fabs(*dest) > 2000))
		  fprintf (stderr,"(%f %f %f) ", *dest, w_1, w_2);*/
	     }
	}
    }

  compute_low_high_hys_thresh (edge_ptr, slope_ptr, LOW_HYSTERESIS, 
			       HIGH_HYSTERESIS, 
			       &LOW_TH, &HIGH_TH);

  for (i=0; i < img_size; i++) {
    if ((*(edge_ptr+i) > 0) && (fabs(*(slope_ptr+i)) < LOW_TH))
      *(edge_ptr+i) = 0;
  }
  /*save_edges(edge_ptr, "temp2.edges");*/
  /*getchar();*/
  /* removing the one pixel edges and filling one pixel gaps.
     and removing the false "Ringing" edges */
  /*fprintf(stderr, " (%f) ", max_val);*/

  for( i= 0; i<img_height ; i++)
    {  for (j=0; j<img_width; j++)
      { 
	dest1 = edge_ptr + j + i * img_width;
	dest =  slope_ptr + j + i * img_width;
	if ((i==0)||(i== img_height-1)||(j==0)||(j==img_width-1))
	  {*dest1 =0;}
	else if (((*dest1 > 0) && (*(dest1 - 1) == 0) && (*(dest1 +1) == 0)
		  && (*(dest1 - img_width)== 0) && (*(dest1 + img_width) == 0)
		  && (*(dest1 - img_width -1)==0) &&
		  (*(dest1 - img_width +1)==0)
		  && (*(dest1 + img_width -1)==0) &&
		  (*(dest1 + img_width +1)==0))
		 ||((*dest1 > 0) && (*(dest1 - 1)> 0)
		    && (*(dest1 - img_width) > 0))
		 ||((*dest1 > 0) && (*(dest1 +1)> 0)
		    && (*(dest1 - img_width)> 0))
		 ||((*dest1 > 0) &&  (*(dest1 - 1)> 0)
		    && (*(dest1 + img_width) > 0))
		 ||((*dest1 > 0) &&  (*(dest1 + 1)> 0)
		    && (*(dest1 + img_width) > 0)))
	  { *dest1 = 0;}
	else if ((*dest1 == 0) &&
		   ((*(dest1 + 1) > 0)&&(*(dest1 - 1) > 0) &&
		    (*(dest1-img_width-1)==0)&&(*(dest1-img_width+1)==0) &&
		    (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		    (*(dest1-img_width)==0)&&(*(dest1+img_width)==0)))
	    {*dest1 = *(dest1 + 1);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1)==0) &&
		   (*(dest1-img_width-1)==0)&&(*(dest1-img_width+1)==0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)>0)&&(*(dest1+img_width)>0))
	    {*dest1 = *(dest1 + img_width);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1)==0) &&
		   (*(dest1-img_width-1)==0)&&(*(dest1-img_width+1)>0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)>0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)==0))
	    {*dest1 = *(dest1 - img_width + 1);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1)==0) &&
		   (*(dest1-img_width-1)>0)&&(*(dest1-img_width+1)==0) &&
		   (*(dest1+img_width+1)>0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)==0))
	    {*dest1 = *(dest1 + img_width + 1);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1) > 0) &&
		   (*(dest1-img_width-1)==0)&&(*(dest1-img_width+1)>0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)==0))
	    {*dest1 = *(dest1 - 1);}
	  else if ((*(dest1 + 1) > 0)&&(*(dest1 - 1)== 0) &&
		   (*(dest1-img_width-1)>0)&&(*(dest1-img_width+1)==0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)==0))
	    {*dest1 = *(dest1 + 1);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1)==0) &&
		   (*(dest1-img_width-1)>0)&&(*(dest1-img_width+1)==0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)>0))
	    {*dest1 = *(dest1 + img_width);}
	  else if ((*(dest1 + 1)==0)&&(*(dest1 - 1)==0) &&
		   (*(dest1-img_width-1)==0)&&(*(dest1-img_width+1)>0) &&
		   (*(dest1+img_width+1)==0)&&(*(dest1+img_width-1)==0) &&
		   (*(dest1-img_width)==0)&&(*(dest1+img_width)>0))
	    {*dest1 = *(dest1 + img_width);}
	
	
      }
    }
  /* Save the edge imgae */
  /*save_edges(edge_ptr, convlfilename);*/
  /*save_edges(edge_ptr, "temp2.edges");*/
  if ((convl_ptr = (float *)calloc((img_size), sizeof(float))) == NULL)
   {
      fprintf(stderr, "Problems allocating memory for convolved image") ;
      exit(1) ;  
    }
  for (i=0; i< img_height;i++)
    {
      for (j=0; j< img_width; j++)
	{ *(convl_ptr + i*img_width + j) = 
	   *(derivx_ptr + (i+gksize)*expnd_width + j + gksize);
	}
    }
  free(derivx_ptr);

  if ((dir_ptr = (float *)calloc((img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for convolved image") ;
      exit(1) ; 
    }
  for (i=0; i< img_height;i++)
    {
      for (j=0; j< img_width; j++)
	*(dir_ptr + i*img_width + j) =
	  *(derivy_ptr + (i+gksize)*expnd_width + j + gksize);
    }

  
  free(derivy_ptr);

}   


/*******************************************************************/
int derivy_convl (UCHAR *imgptr, float *outptr)
{ 
  register int i, j ;

  display_main(" Doing column convolution by the Canny operator.");
  /*fprintf(stderr, "\n OPERATOR : %d", OPERATOR);*/

  b1 = 3 * (exp (- beta * tau_y));
  b2 = -3 * (exp (-2 * beta * tau_y));
  b3 = exp(-3 * beta * tau_y);
  a0 = (alpha + 2)/ (beta * beta);
  a1 = ((tau_y * (alpha + 2))/ beta + 0.5 * tau_y * tau_y * (alpha + 1) - 
	2 * (alpha + 2) / (beta * beta)) * exp (- beta * tau_y);
  a2 = (- (tau_y * (alpha + 2))/ beta + 0.5 * tau_y * tau_y * (alpha + 1) 
	+ (alpha + 2) / (beta * beta)) * exp (-2 * beta * tau_y);
  a01 = (0.5 * beta * tau_y * tau_y * (alpha + 1) + tau_y) 
    * exp( - beta * tau_y);
  a02 = (0.5 * beta * tau_y * tau_y * (alpha + 1) - tau_y) 
    * exp(-2 *  beta * tau_y);
  
  for (j=0; j < expnd_width; j++) {
    for (i=0; i < expnd_height; i++) { 
      *(temp_array+i) = *(imgptr + j + (i * expnd_width)); }
    if(opt_convl(temp_array, expnd_height) == NOOK) {
      fprintf(stderr, "Failure in projection recursive filtering.\n") ;
      exit(1) ; }
    for (i=0; i < expnd_height; i++) {
      *(outptr + j + (i * expnd_width)) = *(temp_array+i);}
  }

  display_main(" Doing row convolution by the Canny  projection function.");
  /* the constant are defined so as to works with different tau_x's and
     tau_y's*/

  b1 = 3 * (exp (- beta * tau_x));
  b2 = -3 * (exp (-2 * beta * tau_x));
  b3 = exp(-3 * beta * tau_x);
  a0 = (alpha + 2)/ (beta * beta);
  a1 = ((tau_x * (alpha + 2))/ beta + 0.5 * tau_x * tau_x * (alpha + 1) - 
	2 * (alpha + 2) / (beta * beta)) * exp (- beta * tau_x);
  a2 = (- (tau_x * (alpha + 2))/ beta + 0.5 * tau_x * tau_x * (alpha + 1) 
	+ (alpha + 2) / (beta * beta)) * exp (-2 * beta * tau_x);
  a01 = (0.5 * beta * tau_x * tau_x * (alpha + 1) + tau_x) 
    * exp( - beta * tau_x);
  a02 = (0.5 * beta * tau_x * tau_x * (alpha + 1) - tau_x) 
    * exp(-2 *  beta * tau_x);
  
    
  for (i=0; i < expnd_height; i++)
    {
      for (j=0; j < expnd_width; j++)
	{ 
	  *(temp_array+j) = (float)(*(outptr + j + (i * expnd_width)));
	}
      if(projec_convl(temp_array, expnd_width) == NOOK)
	{
	  fprintf(stderr, "Failure in projection recursive filtering.\n") ;
	  exit(1) ;  /* return(NOOK) ; */
	}
      for (j=0; j < expnd_width; j++)
	{
	  *(outptr + j + (i * expnd_width)) = *(temp_array+j);  
	}
    }

}



/****************************************************************/
int derivx_convl (UCHAR *imgptr, float *outptr)
{ 
  register int i, j ;
  char file2[100];

  display_main(" Doing row convolution by Canny operator.");
  b1 = 3 * (exp (- beta * tau_x));
  b2 = -3 * (exp (-2 * beta * tau_x));
  b3 = exp(-3 * beta * tau_x);
  a0 = (alpha + 2)/ (beta * beta);
  a1 = ((tau_x * (alpha + 2))/ beta + 0.5 * tau_x * tau_x * (alpha + 1) - 
	2 * (alpha + 2) / (beta * beta)) * exp (- beta * tau_x);
  a2 = (- (tau_x * (alpha + 2))/ beta + 0.5 * tau_x * tau_x * (alpha + 1) 
	+ (alpha + 2) / (beta * beta)) * exp (-2 * beta * tau_x);
  a01 = (0.5 * beta * tau_x * tau_x * (alpha + 1) + tau_x) 
         * exp( - beta * tau_x);
  a02 = (0.5 * beta * tau_x * tau_x * (alpha + 1) - tau_x) 
         * exp(-2 *  beta * tau_x);
  a0 *= beta*beta;
  a1 *= beta*beta;
  a2 *= beta*beta;

  for (i=0; i < expnd_height; i++)
    {
      for (j=0; j < expnd_width; j++)
	{ 
	  *(temp_array+j) = (float)(*(imgptr + j + (i * expnd_width)));
	}

      if(opt_convl(temp_array, expnd_width) == NOOK)
	{
	  fprintf(stderr, "Failure in projection recursive filtering.\n") ;
	  exit(1) ;  /* return(NOOK) ; */
	}
  

      for (j=0; j < expnd_width; j++)
	{
	  *(outptr + j + (i * expnd_width)) = *(temp_array+j);
	}
    }

  
  display_main(" Doing column convolution by Canny projection function.");
  b1 = 3 * (exp (- beta * tau_y));
  b2 = -3 * (exp (-2 * beta * tau_y));
  b3 = exp(-3 * beta * tau_y);
  a0 = (alpha + 2)/ (beta * beta);
  a1 = ((tau_y * (alpha + 2))/ beta + 0.5 * tau_y * tau_y * (alpha + 1) - 
	2 * (alpha + 2) / (beta * beta)) * exp (- beta * tau_y);
  a2 = (- (tau_y * (alpha + 2))/ beta + 0.5 * tau_y * tau_y * (alpha + 1) 
	+ (alpha + 2) / (beta * beta)) * exp (-2 * beta * tau_y);
  a0 *= beta*beta;
  a1 *= beta*beta;
  a2 *= beta*beta;
  a01 = (0.5 * beta * tau_y * tau_y * (alpha + 1) + tau_y) 
    * exp( - beta * tau_y);
  a02 = (0.5 * beta * tau_y * tau_y * (alpha + 1) - tau_y) 
    * exp(-2 *  beta * tau_y);
  
  for (j=0; j < expnd_width; j++)
    {
      for (i=0; i < expnd_height; i++)
	{ 
	  *(temp_array+i) = *(outptr + j + (i * expnd_width));
	}
      if(projec_convl(temp_array, expnd_height) == NOOK)
	{
	  fprintf(stderr, "Failure in projection recursive filtering.\n") ;
	  exit(1) ;  /* return(NOOK) ; */
	}
      for (i=0; i < expnd_height; i++)
	{
	  *(outptr + j + (i * expnd_width)) = *(temp_array+i);  
	}
    }
}

/****************************************************************/
int projec_convl (float *vector,  int size)
{  
  register int i;
     
  if (OPERATOR == 2) {
    canny_proj_convl(vector, size);
    return;
  }
  for (i=0; i < size; i++)
    { *(temp_array1 +i)= *(vector + i);
    }
  for (i=0; i < size; i++)
    { *(temp_array2 + (size - i - 1))= *(vector + i);
    }
  
  if(half_projec_convl(temp_array1, size) == NOOK)
    {
      fprintf(stderr, 
	      "Failure in projection recursive filtering of a vector(1).\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  
  if(half_projec_convl(temp_array2, size) == NOOK)
    {
      fprintf(stderr, 
	      "Failure in projection recursive filtering of a vector(2).\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  
  for (i=0; i < size; i++)
    { *(vector + i) = *(temp_array1 + i) + *(temp_array2 + size - i - 1) 
	- (a0 * *(vector + i));
    }
}

/****************************************************************************/
int opt_convl (float *vector,  int size)
{  
  register int i;
  
  if (OPERATOR == 2) {
    canny_opt_convl(vector, size);
    return;
  }  
  for (i=0; i < size; i++)
    { *(temp_array1 +i)= *(vector + i);
    }
  for (i=0; i < size; i++)
    { *(temp_array2 + (size - i - 1))= *(vector + i);
    }
  
  if(half_opt_convl(temp_array1, size) == NOOK)
    {
      fprintf(stderr, 
	      "Failure in projection recursive filtering of a vector(1).\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  
  if(half_opt_convl(temp_array2, size) == NOOK)
    {
      fprintf(stderr, 
	      "Failure in projection recursive filtering of a vector(2).\n") ;
      exit(1) ;  /* return(NOOK) ; */
    }
  for (i=0; i < size; i++)
    { *(vector + i) = *(temp_array1 + i) - *(temp_array2 + size - i - 1);
    }
}
/********************************************************************/
float output[5000];

int canny_opt_convl (float *vector,  int size)
{  
  register int i, j;
  
  if (size > 5000) {
    DEBUG
      fprintf(stderr, "\n Image size is too large");
    exit;
  }
  
  
  for (i=0; i < size; i++) {
    output[i] = 0.0;
    for (j=-W; j <=W; j++) {
      if (((j+i) >= 0) && ((j+i) < size)) {
	output[i] += canny_op[j+W] * vector[j+i];
      }
    }
  }
  for (i=0; i < size; i++) {
    vector[i] = output[i];
  }
}


/********************************************************************/

int canny_proj_convl (float *vector,  int size)
{  
  register int i, j;
  
  if (size > 5000) {
    DEBUG
      fprintf(stderr, "\n Image size is too large");
    exit;
  }
  
  
  for (i=0; i < size; i++) {
    output[i] = 0.0;
    for (j=-W; j <=W; j++) {
      if (((j+i) >= 0) && ((j+i) < size)) {
	output[i] += canny_sm[j+W] * vector[j+i];
      }
    }
  }
  for (i=0; i < size; i++) {
    vector[i] = output[i];
  }
}


/********************************************************************/


int canny_mask (float sigma, int mrc) 
{
  float a, w, b, sum = 0.0;
  float a1, a2, a3, a4, dt, t;
  int i, j;

  float coeff[7][5] = {{24.59550, 0.12250, 63.97566, 0.15, 4.21},
		       {12.47120, 0.38284, 31.26860, 0.3,  2.87},
		       {7.85869,  2.62856, 18.28800, 0.5,  2.13},
		       {5.06500,  2.56770, 11.06100, 0.8,  1.57},
		       {3.45580,  0.07161, 4.80684,  1.0,  1.33},
		       {2.05220,  1.56939,  2.91540, 1.2,  1.12},
		       {0.00297,  3.50350,  7.47700, 1.4,  0.75}};

  float A[7][4] =    {{4.1447324574546559e-09,  -5.3125372670207223e-10,   3.2147069408673610e+02,  -9.9999999946874629e-01},
		      {1.1448154283245602e-04,  -5.0237960566262104e-05,   4.9103043640146709e+01,  -9.9994976203943386e-01},
		      {-1.1988842971055013e-03,  -2.3143241098873575e-04,   3.9702722595093252e+00,  -9.9976856758901123e-01},
		      {-1.4260822017760764e-02,  -1.6003357995773733e-03,   2.3557381759336069e+00,  -9.9839966420042259e-01},
		      {1.6213161968627832e+00,  -1.4935950384115138e-01,   3.1661090544493206e+01,  -8.5064049615884851e-01},
		      {-1.4867687163692719e-01,  -2.0875534760151188e-01,  1.2446539383288555e+00,  -7.9124465239848818e-01},
		      {-2.6160798806796680e+02,   3.4893306749433236e+02,   2.6314969424153134e+02,  -3.4993306749433259e+02}};


  

  if (mrc > 6) mrc = 6;
  a=coeff [mrc][0];
  w=coeff [mrc][1];
  b=coeff [mrc][2]; 

  dt = 1/(4*sigma); W = (int) (1/dt);
  if (W < 2) {W = 4; dt = 0.25;}
  
  if (W > 500) {
    fprintf(stderr, "\n The Canny operator mask is too large");
    DEBUG
      exit;
  }
  
  /*a1 = (a*(b-a)*sin(2*w) - a*w*cos(2*w) + (-2*w*w*sinh(a) 
    +2*a*a*exp(-a))*sin(w) + 2*a*w*sinh(a)*cos(w)
    +w*exp(-2*a)*(a+b) - b*w)/(4*(w*w*sinh(a)*sinh(a) - a*a*sin(w)*sin(w)));
    
    a2 = (a*(b-a)*cos(2*w) + a*w*sin(2*w) -2*a*w*cosh(a)
    *sin(w) - 2*w*w*sinh(a)*cos(w) + 2*w*w*exp(-a)*sinh(a)
    + a*(a-b))/(4*(w*w*sinh(a)*sinh(a) - a*a*sin(w)*sin(w)));
    
    a3 = (-a*(b+a)*sin(2*w) + a*w*cos(2*w) + (2*w*w*sinh(a)
    + 2*a*a*exp(a))*sin(w) + 2*a*w*sinh(a)*cos(w)
    + w*exp(2*a)*(b-a) - b*w)/(4*(w*w*sinh(a)*sinh(a)-a*a*sin(w)*sin(w)));
    
    a4 = (-a*(b+a)*cos(2*w) - a*w*sin(2*w) + 2*a*w*cosh(a)
    *sin(w) + 2*w*w*sinh(a)*cos(w) - 2*w*w*exp(a)*sinh(a)
    + a*(a-b))/(4*(w*w*sinh(a)*sinh(a)-a*a*sin(w)*sin(w)));
    */
  a1 = A[mrc][0]; a2 = A[mrc][1]; a3 = A[mrc][2]; a4 = A[mrc][3];

   t = 1;
   /*fprintf(stderr, "\n %f %f %f %f %f", a1, a2, a3, a4, exp(a*t)*
     (a1*sin(w*t) + a2*cos(w*t))+
     exp(-a*t)*(a3*sin(w*t) + a4*cos(w*t)) + 1);*/

   for (t=0, i=W; t <=1; t = t + dt, i++) {
     canny_op[i] = exp(a*t)*(a1*sin(w*t) + a2*cos(w*t)) + 
		   exp(-a*t)*(a3*sin(w*t) + a4*cos(w*t)) + 1;
   }

  for (i=0, j=W; i >= -W; i--, j++) {
    sum += canny_op[i];
    canny_op[i+W] = - canny_op[j];
  }
  
  for (i=-W; i<= W; i++) {
    if (i== -W) canny_sm[0] = canny_op[0];
    /*else  if (i== W) canny_sm[0] = 0.0;*/
    else canny_sm[i+W] = canny_sm[i+W-1] + canny_op[i+W];
    
    /*fprintf(stderr, "\n %f %f", canny_op[i+W], canny_sm[i+W]);*/
  }
}

/*********************************************************************/
/**This is the function to compute the causal filtering by the 
   projection function */

int half_projec_convl  (float *vector, int size)
{
  int i;
  float x_1, x_2, temp;
  
  x_1 = *(vector + 1);
  x_2 = *(vector + 2);
  *vector = 0.0;
  *(vector + 1) = 0.0;
  *(vector + 2) = 0.0;
  for (i=3; i < size; i++)
    {temp = ( a0 * *(vector + i)) + (a1 * x_1) + (a2 * x_2) 
       + (b1 * *(vector +i - 1)) + (b2 * *(vector + i - 2)) 
	 + (b3 * *(vector + i -3));
     x_2 = x_1;
     x_1 = *(vector + i);
     *(vector + i) = temp;
   }
}

/*********************************************************************/
/**This function computes the causal operator filtering **/


int half_opt_convl  (float *vector,  int size)
{
  int i;
  float x_1, x_2, temp;
  
  x_1 = *(vector + 1);
  x_2 = *(vector + 2);
  *vector = 0.0;
  *(vector + 1) = 0.0;
  *(vector + 2) = 0.0;
  for (i=3; i < size; i++)
    {temp = (a01 * x_1) + (a02 * x_2) 
       + (b1 * *(vector +i - 1)) 
	 + (b2 * *(vector + i - 2)) + (b3 * *(vector + i -3));
     x_2 = x_1;
     x_1 = *(vector + i);
     *(vector + i) = temp;
     
   }
}



