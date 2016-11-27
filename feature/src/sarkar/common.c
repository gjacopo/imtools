/* This file contains common functions that are used by other files */
  

#include "header.h"
/******************************************************************/
float my_pow (float x, float d)
{
  return(pow(fabs(x), d));
}

/******************************************************************/
int slstore3(struct scale_data *i, struct scale_data **last, 
	     struct scale_data **start)
{
  if (!*last) /* first item in the list */
    {
      *last = i;
      *start = i;
    } 
  else (*last)->next = i;
  i->next = NULL;
  *last = i;
}

/******************************************************************/

float cot (x)
     float x;
{
  float y;
  y = 1.0 / tan(x);
  return (y);
}

/***********************************************************************/
     
int save_edges(UCHAR *zcptr, char  filename[])
{
  FILE *outfile1;
  int i, imax;
  char file1[100];

  /*fprintf(stderr, "\n Saving edge image in %s...\n", filename) ;*/
  strcpy(file1, filename) ;

  if ((outfile1 = fopen(file1, "w")) == NULL)
    {
      fprintf(stderr, "Cannot open file, %s for writing.\n", file1) ;
      return(NOOK) ;
    }
  fprintf(outfile1, "P5\n%d %d\n255\n", img_width, img_height);

  imax = img_width * img_height ;
  for (i = 0; i < img_width; i++)
    *(zcptr + i) = 255;
  for ( i = imax; i; i--, zcptr++ )
    putc(*zcptr, outfile1) ;  
  fclose(outfile1) ;
  return(OK) ;
}

/************************************************************************/
    
int save_convl(float *convlptr, char  file2[])
{
  FILE  *outfile2 ;
  float *convl, convert, max_value ;
  int i, j;
  UCHAR newvalue ;



  if ((outfile2 = fopen(file2, "w")) == NULL)
    {
      fprintf(stderr, "Cannot open file, %s for writing.\n", file2) ;
      return(NOOK) ;
    }
  
  convl = convlptr ;

  fprintf(outfile2, "P5\n%d %d\n255\n", img_width, img_height);
  max_value = 0.0 ;
  for ( i = 0; i< img_height; i++)
    {
      for (j=0; j < img_width; j++)
	{
	  if ( fabs( *(convl + j +gksize + ((i+gksize) * expnd_width))) 
	      > max_value )
	    max_value = fabs( *(convl + j +gksize + 
				((i+gksize) * expnd_width))) ;
	  /* find the max value of the image */ 
	}
    }
  fprintf(stderr, "max_value= %f\n", max_value) ;
  convert = 255.0 / max_value ;
  convl = convlptr ;
  /*fprintf(outfile2, "%f\n", max_value) ;*/
  for ( i = 0; i< img_height; i++)
    {
      for (j=0; j < img_width; j++)
	{
	  if ( *(convl + j +gksize + ((i+ gksize) * expnd_width)) < 0.0)
	    newvalue =(UCHAR)(int)( - *(convl + j +gksize 
					+ ((i+gksize) * expnd_width)) 
				   * convert ) ;
	  else
	    newvalue =(UCHAR)(int)( *(convl + j +gksize 
				      + ((i+gksize) * expnd_width)) 
				   * convert ) ;

	  putc(newvalue, outfile2);
	}
    }
  fclose(outfile2) ;
  return(OK) ;
}


/**************************************************************************/
float angle_diff (float theta1, float theta2)
/* anti clockwise => positive difference */
{
  float temp;
  temp = theta1 - theta2;
  if (fabs(temp) < 3.142)
    return(temp);
  else
    {if (temp < 0)
       return(temp + 6.284);
     else
       return(temp - 6.284);
   }
}

/************************************************************************
  Code to maintain an integer linked list ******/

int length (struct int_list *list)

{
  int count;
  struct int_list *temp;

  count = 0;
  temp = list;
  while (temp)
    {count++;
     temp = temp->next;
   }
  return (count);
}

/*************************************************************************/
int delete_list (struct int_list *head)
{
  if (head != NULL) {
    delete_list ( head->next);
    free (head);
  }
}

/*************************************************************************/
int append_element (struct int_list **list, int value)
{
  struct int_list *temp, *temp1;
  
  if ((*list) == NULL) { /* first element */
    temp = (struct int_list *) malloc (sizeof(Int_list));
    temp->value = value;
    temp->next = NULL;
    (*list) = temp;
  }
  else {
    temp = (*list);
    while (temp->next) {
      temp = temp->next;
    }
    temp1 = (struct int_list *) malloc (sizeof(Int_list));
    temp1->value = value;
    temp1->next = NULL;
    temp->next = temp1;
  }
  
}

/*************************************************************************/

int print_list (struct int_list *list)
{
  if (list != NULL) {
    fprintf(stderr, " %d", list->value);
    print_list (list->next);
  }
}

/*************************************************************************/
int cons_element (struct int_list **list, int value)
{
  struct int_list *temp;
  temp = (struct int_list *) malloc (sizeof(Int_list));
  temp->value = value;
  temp->next = *list;
  *list = temp;
}
/**************************************************************************/

int max(int x, int y)
{
  if (x > y)
    {return (x);}
  else
    {return (y);}
}

/******************************************************************/

float fmax(double x, double y)
{
  if (x > y)
    {return (x);}
  else
    {return (y);}
}

/******************************************************************/
float fmin(double x, double y)
{
  if (x <= y)
    {return (x);}
  else
    {return (y);}
}

/******************************************************************/
  
int min(int x, int y)
{
  if (x <= y)
    {return (x);}
  else
    {return (y);}
}

 
/******************************************************************/

int consistent (int dir1, int dir2)
{
  if (((dir1 == 0) && ((dir2 == 7) || (dir2 == 0) || (dir2 == 1)))||
      ((dir1 == 1) && ((dir2 == 0) || (dir2 == 1) || (dir2 == 2)))||
      ((dir1 == 2) && ((dir2 == 1) || (dir2 == 2) || (dir2 == 3)))||
      ((dir1 == 3) && ((dir2 == 2) || (dir2 == 3) || (dir2 == 4)))||
      ((dir1 == 4) && ((dir2 == 3) || (dir2 == 4) || (dir2 == 5)))||
      ((dir1 == 5) && ((dir2 == 4) || (dir2 == 5) || (dir2 == 6)))||
      ((dir1 == 6) && ((dir2 == 5) || (dir2 == 6) || (dir2 == 7)))||
      ((dir1 == 7) && ((dir2 == 6) || (dir2 == 7) || (dir2 == 0))))
    return(1);
  else
    return(0);
}
 
/******************************************************************/

int reverse_list (struct c_element **start)
{
  struct c_element *last, *temp, *temp1, *temp2;
  int length;
  
  temp = *start;
  temp1 = temp->next;
  while(temp1){
    temp2 = temp1->next;
    temp1->next = temp;
    temp = temp1;
    temp1 = temp2;
  }
  (*start)->next = NULL;
  *start = temp;
}



