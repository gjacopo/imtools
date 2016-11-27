/* This file contains code to do thresholding on the edges.
   It accepts a chain code of edge points, a high threshold. It keeps 
   edge chains with aleast 20% points above the high threshold. */

#include "header.h"

/*********************************************************************************/
int threshold_edges (Chain *edges, float thresh)
{
  Chain *temp;
  Element *temp1;
  int count, len;
  
  temp = edges;
  while (temp) {
    if (temp->delete_p != 2) {
      /* do not consider edge chains deleted by edge focusing */
      temp1 = temp->start; count = len = 0;
      while (temp1) {
	if (fabs(temp1->slope) > thresh) count++;
	len++;
	temp1 = temp1->next;
      }
      if (count < 0.2*len) temp->delete_p = 1;
      else temp->delete_p = 0;
    }
    temp = temp->next;
  }
}

/*********************************************************************************/

int threshold_edges_length (Chain *edges, int length_thresh)
{
  Chain *temp;
  Element *temp1;
  int count, len;
  
  temp = edges;
  while (temp) {
    if (temp->delete_p == 0) {
      /* do not consider edge chains deleted by edge focusing */
      temp1 = temp->start; count = len = 0;
      while (temp1) {
	len++;
	temp1 = temp1->next;
      }
      if (len < length_thresh) temp->delete_p = 1;
      else temp->delete_p = 0;
    }
    temp = temp->next;
  }
}

/*********************************************************************************/
int floatcompare(float *i, float *j)
{
  if (*i > *j)
    return (1);
  if (*i < *j)
    return (-1);
  return (0);
}
/****************************************************************************/
static float *a;
static int n;

/****************************************************************************/

int compute_low_high_hys_thresh (UCHAR *edges, float *slopes, float LH, 
				 float HH, float *LTH, float *HTH)
{
  int i;
  
  n = 0;
  if ((a = (float *)calloc((img_size), sizeof(float))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for computing hysteresis thresh.\n");
      exit(1) ;  
    }

  for (i=0; i < img_size; i++) {
    if (*(edges+i) > 0)   a[n++] = *(slopes + i);
  }

  qsort((void *) a, n, sizeof(float), floatcompare);

  *LTH = a[(int) (LH*(n-1))];   *HTH = a[(int) (HH*(n-1))];
  
  /*fprintf(stderr, "\n Low Thresh: %f High Thresh: %f", *LTH, *HTH);*/
}
/****************************************************************************/

int recompute_low_high_hys_thresh (float LH, float HH, float *LTH, float *HTH)
{

  *LTH = a[(int) (LH*(n-1))];   *HTH = a[(int) (HH*(n-1))];
  
  /*fprintf(stderr, "\n Low Thresh: %f High Thresh: %f", *LTH, *HTH);*/
}

/***************************************************************************/
