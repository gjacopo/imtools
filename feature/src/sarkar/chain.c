/* Has code to create the chain code the edge points */

/**************************************************************************/
extern ATTS_FILE_P;
#include "header.h"

void initialize_element (Element *temp)
{
  temp->curvature = 0.0;
  temp->dir = 0.0;
  temp->tan_dir = 0.0;
  temp->grad_dir = 0.0;
  temp->slope_ratio = 0.0;
  temp->x = 0;
  temp->y = 0;
  temp->slope = 0.0;
  temp->mag_plus = 0.0;
  temp->width_plus = 0.0;
  temp->mag_minus = 0.0;
  temp->width_minus = 0.0;
}

int write_chain (char filename[100], struct chain **start)
{
  int  x, y,flag, flag1,id;
  struct chain *temp;
  struct c_element *temp1;
  register int i, count;
  char file1[100];
  UCHAR *edge_ptr;
 
  strcpy(file1, filename) ;
  strcat(file1, ".edge") ;
    
  if ((edge_ptr = (UCHAR *)calloc((UINT)(img_size), sizeof(UCHAR))) == NULL)
    {
      fprintf(stderr, "Problems allocating memory for edge array ") ;
      exit(1) ;  
    }
  for (i=0;i<img_size;i++)
    {*(edge_ptr + i) = 255;}
  temp = *start;
  count = 0;
  while (temp) {
    /*fprintf(stderr, " %d", count++);
      fprintf(stderr, ".");*/
    if (temp->delete_p == 0) {
      temp1 = temp->start;
      x = temp->start_x;
      y = temp->start_y;
      temp1 = temp1->next;
      while (temp1){
	x = temp1->x;
	y = temp1->y;
	if (ATTS_FILE_P == 1) 
	  *(edge_ptr + x + y*img_width) = 0; /*160+temp1->tan_dir;*/
	else 
	  *(edge_ptr + x + y*img_width) = 0;

	temp1 = temp1->next;
      }
    }
    temp = temp->next;
  }
  temp = *start;
  while (temp) {
    if (temp->delete_p == 0) {    
      *(edge_ptr + temp->start_x + (temp->start_y)*img_width) = 0;
      if (ATTS_FILE_P == 1)
	*(edge_ptr + temp->end_x + (temp->end_y)*img_width) = 0; /*50;*/
      else 
	*(edge_ptr + temp->end_x + (temp->end_y)*img_width) = 0;
    }
    temp = temp->next;
  }
  save_edges(edge_ptr, file1);
  free(edge_ptr);
}


/************************************************************************/
void slstore(struct chain *i, struct chain **last, struct chain **start)
{
  i->delete_p = 0;
  if (!*last) /* first item in the list */
    {
      *last = i;
      *start = i;
    } 
  else (*last)->next = i;
  i->next = NULL;
  *last = i;
}
/************************************************************************/
void insert_Element (Element **head, Element *ei)
{
  Element *temp;

  if (!*head) *head = ei;
  else {
    temp = *head;
    while (temp->next) temp = temp->next;
    temp->next = ei;
  }
}
/************************************************************************/
int num_of_neighbor (UCHAR *point)
{
  int temp;
  
  temp = ( *(point  + 1) + *(point  - 1) + *(point + img_width)
	  + *(point  - img_width) + *(point  - img_width + 1)
	  + *(point  + img_width + 1) + *(point  + img_width - 1)
	  + *(point  - img_width - 1));
  temp = (int)(temp/240);
  return (temp);
}
/***********************************************************************/
struct chain *chain_coder ()
{
  int i, j;
  struct chain *temp;
  struct c_element *temp1, *temp2;
  UCHAR *pixel, *nextp;
  struct chain *start;
  struct chain *last;
  int x, y, count, ii, jj;

  start = NULL;
  last = NULL;
  display_main("Chain coding the edges...");
  count = 0;
  pixel = edge_ptr;
  for (i=0;i<img_height;i++) {
    for (j=0; j<img_width;j++, pixel++) {
      
      if ((*pixel > 0) && (num_of_neighbor(pixel)== 1)) {
	temp = (struct chain *)malloc(sizeof(Chain));
	   if (!temp) 
	     {fprintf( stderr, "Out of memory \n");
	      exit (1);
	      }
	   /*fprintf(stderr, "%d ", count);*/
	   temp->id = count++;
	   temp->start_x = j;
	   temp->start_y = i;
	   temp1 = (struct c_element *)malloc(sizeof(Element));
	   if (!temp1) 
	     {fprintf( stderr, "Out of memory \n");
	      exit (1);
	      }
	   initialize_element (temp1);
	   temp->start = temp1;
	   temp1->dir  = *pixel - 240;
	   temp1->x = j; temp1->y = i;
	   get_attbs (pixel, temp1);
	   grow_seed(pixel, &temp1, j, i);
	   while (temp1->next) temp1 = temp1->next;
	   temp->end_x = temp1->x;
	   temp->end_y = temp1->y;
	   slstore( temp, &last, &start);
	 }
      }
    }
  /* This part handles rest of the the branch points or joined closed
     contours */
  pixel = edge_ptr;
  for (i=0;i<img_height;i++)
    {for (j=0; j<img_width;j++, pixel++)
       {if ((*pixel > 0) && (num_of_neighbor(pixel) > 2))
	 {*pixel = 0; 
	  for (ii= -1;ii < 2;ii++)
	    {for (jj= -1;jj < 2;jj++)
	       {if ((*(pixel + ii * img_width + jj) > 0)
		    && (num_of_neighbor(pixel+ ii * img_width + jj) ==1)){
		 temp = (struct chain *)malloc(sizeof(Chain));
		 if (!temp) 
		   {fprintf( stderr, "Out of memory \n");
		    exit (1);
		  }
		 /*fprintf(stderr, "%d ", count);*/
		 temp->id = count++;
		 temp->start_x = j + jj;
		 temp->start_y = i + ii;
		 
		 temp1 = (struct c_element *)malloc(sizeof(Element));
		 if (!temp1) 
		   {fprintf( stderr, "Out of memory \n");
		    exit (1);
		  }
		 initialize_element (temp1);
		 temp->start = temp1;
		 temp1->dir  = *(pixel+ ii * img_width + jj) - 240;
		 temp1->x = j + jj; temp1->y = i + ii;
		 get_attbs (pixel, temp1);

		 grow_seed((pixel+ ii * img_width + jj), &temp1, j + jj, i+ii);
		 while (temp1->next) temp1 = temp1->next;
		 temp->end_x = temp1->x;
		 temp->end_y = temp1->y;
		 slstore( temp, &last, &start);
	       }
	      }
	    }
	 }
      }
   }
  /* This part lists al the remaining closed contours */
  display_main(" Listing closed contours ");
  pixel =edge_ptr;
  for (i=0;i<img_height;i++)
    {for (j=0; j<img_width;j++, pixel++)
       {if ((*pixel > 0) && (num_of_neighbor(pixel)== 2))
	  {temp = (struct chain *)malloc(sizeof(Chain));
	   if (!temp) 
	     {fprintf( stderr, "Out of memory \n");
	      exit (1);
	      }
	   /*fprintf(stderr, "%d ", count);*/
	   temp->id = count++;
	   temp->start_x = j;
	   temp->start_y = i;
	   temp1 = (struct c_element *)malloc(sizeof(Element));
	   if (!temp1) 
	     {fprintf( stderr, "Out of memory \n");
	      exit (1);
	      }
	   initialize_element (temp1);
	   temp2 = (struct c_element *)malloc(sizeof(Element));
	   if (!temp2) 
	     {fprintf( stderr, "Out of memory \n");
	      exit (1);
	      }
	   temp->start = temp1;
	   initialize_element (temp2);
	   temp1->dir = *pixel - 240;
	   temp1->x = j;
	   temp1->y = i;
	   get_attbs (pixel, temp1);

	   x = j; y = i;
	   *pixel = 0;
	   if (*(pixel + 1) >0) 
	     {nextp = pixel + 1;
	      x++;
	      temp2->dir = 2;}
	   else
	     if (*(pixel - 1) >0) 
	       {nextp = pixel - 1;
		x--;
		temp2->dir= 6;}
	     else
	       if (*(pixel + img_width + 1)) 
		 {nextp = pixel + img_width + 1;
		  x++; y++;
		  temp2->dir = 3;}
	       else
		 if (*(pixel + img_width - 1)) 
		   {nextp = pixel + img_width - 1;
		    x--; y++;
		    temp2->dir = 5;}
		 else
		   if (*(pixel + img_width)) 
			 {nextp = pixel + img_width;
			  y++;
			  temp2->dir = 4;}
	   
	   temp2->x = x;
	   temp2->y = y;
	   get_attbs(nextp, temp2);
	   temp1->next = temp2;
	   temp2->next = NULL;
	   grow_seed(nextp, &temp2, x, y);
	   while (temp2->next) temp2 = temp2->next;
	   temp->end_x = temp2->x;
	   temp->end_y = temp2->y;
	   slstore( temp, &last, &start);
	 }
      }
   }

  free(edge_ptr);

  free(slope_ptr);

  free(convl_ptr);

  free(dir_ptr);


  return(start);
}


/***********************************************************************/
int grow_seed (UCHAR *present,  struct c_element **last, int x, int y)
     
{
  struct c_element *temp;
  UCHAR *nextp;
  int n_nhb;
  int dirn;
  
  n_nhb = num_of_neighbor (present);
  switch (n_nhb)
    {case 1: 
       break;
     default: 
       return;
     }
  *present = 0;
  temp = (struct c_element *)malloc(sizeof(Element));
  initialize_element (temp);
  if (!temp) 
    {fprintf( stderr, "Out of memory \n");
     exit (1);
     }
  if (*(present + 1) >0) 
    {nextp = present + 1;
     x++;
     dirn = 2;}
  else
    if (*(present - 1) >0) 
      {nextp = present - 1;
       x--;
       dirn = 6;}
    else
      if (*(present + img_width + 1)) 
	{nextp = present + img_width + 1;
	 x++; y++;
	 dirn = 3;}
      else
	if (*(present + img_width - 1)) 
	  {nextp = present + img_width - 1;
	   x--; y++;
	   dirn = 5;}
	else
	  if (*(present - img_width + 1)) 
	    {nextp = present - img_width + 1;
	     x++; y--;
	     dirn = 1;}
	  else  
	    if (*(present - img_width - 1)) 
	      {nextp = present - img_width - 1;
	       x--; y--;
	       dirn = 7;}
	    else
	      if (*(present + img_width)) 
		{nextp = present + img_width;
		 y++;
		 dirn = 4;}
	      else
		if (*(present - img_width)) 
		  {nextp = present - img_width;
		   y--;
		   dirn = 0;}
    
  temp->dir = dirn;
  temp->x = x;
  temp->y = y;
  get_attbs(nextp, temp);
  (*last)->next = temp;
  temp->next = NULL;
  grow_seed(nextp, &temp, x, y);

}
  


/************************************************************************/
/* This function calculates the attributes of the edge points like
   slope, mag+, mag-, width+ , width- */

int get_attbs(UCHAR *edge, struct c_element *f_list)
{
  int flag, offset, row, col;
  register int i, j,k;
  float val1, val2;
  float rowinc, colinc;
  int r,c;
  int width;

  
  offset = edge - edge_ptr;
  row = (int)(offset / img_width);
  col = offset - row*img_width;
  
  
  f_list->tan_dir = (*edge) - 240;
  f_list->grad_dir = *(dir_ptr + offset);
  f_list->slope = fabs(*(slope_ptr + offset));
  
  rowinc = sin(f_list->grad_dir); /** Because of inverse co-ord sys*/
  colinc = cos(f_list->grad_dir);
  r = row ; c = col;
  val2 = *(convl_ptr + r*img_width + c);
  f_list->mag_plus = 0.0;
  f_list->mag_minus = 0.0;
  f_list->width_minus = 0.0;
  f_list->width_plus = 0.0;
  for (width=0,flag = 0;
       ((flag==0)&&(r < img_height)&&(c < img_width)&&(r >= 0) &&(c>=0));
       width++)
    { r = (int)(row + width*rowinc);
      c = (int)(col + width*colinc);
      val1 = *(convl_ptr + r*img_width + c);
      /*fprintf(stderr," row= %d col= %d  %f test1.5\n", r, c, val1);*/
      if (val1 < val2)
       {f_list->mag_plus = val2;
	flag = 1;
	if ((rowinc ==0)||(colinc==0))
	  {f_list->width_plus =  (float)width;}
	else 
	  {f_list->width_plus = (sqrt(2.0) * (float)width);}
      }
      val2 = val1;
    }
  r = row ; c = col;
  val2 = *(convl_ptr + r*img_width + c);
  for (width= 0,flag = 0;
       ((flag==0)&&(r < img_height)&&(c< img_width)&&(r >= 0) &&(c>=0));
       width++)
    {r = (int)(row - width*rowinc); 
     c = (int)(col - width*colinc);
     val1 = *(convl_ptr + r*img_width + c);
     if (val1 > val2)
       {f_list->mag_minus = val2;
	flag = 1;
	if ((rowinc ==0)||(colinc==0))
	  {f_list->width_minus = (float)width;}
	else 
	  {f_list->width_minus =( sqrt(2.0)* (float)width);}
      }
     val2 = val1;
   }
  
}

/*********************************************************************/
int chain_length (Element *list)
{
  int count;
  Element *temp;

  count = 0;
  temp = list;
  while (temp)
    {
      count++;
      temp = temp->next;
    }
  return (count);
}


/*********************************************************************/
int chain_count (Chain *list)
{
  int count;
  Chain *temp;

  count = 0;
  temp = list;
  while (temp)
    {
      if (temp->delete_p == 0) count++;
      temp = temp->next;
    }
  return (count);
}

/*********************************************************************/
Chain  *nth_chain (int n, Chain *list)
{
  int count;

  count = 0;
  while (list)
    {
      if (list->delete_p == 0) 
	if (count++ == n) return (list);
      list = list->next;
    }
  return(NULL);
}

/******************************************************************/
void print_chain (Element *list)
{
  Element *temp;

  temp = list;
  while (temp){
    fprintf(stderr, "(%d %d)", temp->x, temp->y);
    temp = temp->next;
  }
}

/*********************************************************/
void element_append (Element **chain, int i, int j)
{
  Element *temp, *temp1;
  
  temp = (Element *) malloc (sizeof (Element));
  initialize_element (temp);
  temp->x = i; temp->y = j; temp->next = NULL;
  
  if (!(*chain)){
    *chain = temp;
  }
  else {
    temp1 = *chain;
    while (temp1->next) temp1 = temp1->next;
    temp1->next = temp;
  }
}

/*************************************************************************/
int delete_Elements (Element *head)
{
  if (head != NULL) {
    delete_Elements ( head->next);
    free (head);
  }
}

/*********************************************************/
int nth_x (int index, Element *list)
{
  int i;
  Element *temp;
  
  temp = list;
  for (i=0; i< index; i++, temp = temp->next) ;
  if (temp) {
    return (temp->x);
  }
  else
    return (-1);   
}

/*********************************************************/
int nth_y (int index, Element *list)
{
  int i;
  Element *temp;
  
  temp = list;
  for (i=0; i< index; i++, temp = temp->next) ;
  if (temp) {
    return (temp->y);
  }
  else
    return (-1);   
}

/*********************************************************/

int add_atts_lists (struct chain *chain_code)
{
  struct chain *temp;
  struct c_element *temp1;
  int length, i, dir, temp_int, id, count;
  float slope, mag_plus, width_plus, mag_minus, width_minus;
  float curvature, grad_dir, slope_0, slope_n, p1, p2;
  double sx, sy;
  char file[100];
  FILE *outfile;

  strcpy(file, convlfilename) ;
  strcat(file, ".atts") ;

  if ((outfile = fopen(file, "w")) == NULL)
    {
      fprintf(stderr, "Cannot open file, %s for writing.\n", file) ;
      return(NOOK) ;
    }
  display_main("  Saving chain code");
  
  fprintf(outfile,"(");
  id = 0;
  temp = chain_code;
  while(temp) {
    if (temp->delete_p == 0) {
      temp1 = temp->start;
      /*fprintf(stderr, " %d", temp->id);*/
      length = 0; slope = 0.0; 
      mag_plus = 0.0; width_plus = 0.0; mag_minus = 0.0;
      width_minus = 0.0; curvature = 0.0; dir = 0; sx = sy = 0;
      while(temp1) {
	sx += temp1->x; sy += temp1->y;
	length++;
	slope = slope + temp1->slope;
	mag_plus = mag_plus + temp1->mag_plus;
	width_plus = width_plus + temp1->width_plus;
	mag_minus = mag_minus + temp1->mag_minus;
	width_minus = width_minus + temp1->width_minus;
	if (fabs(temp1->curvature) < 5.0) /* if curvature has been calcualted*/ 
	  {curvature = curvature + temp1->curvature;}
	if (consistent(temp1->dir, temp1->tan_dir) == 0) dir++;
	temp1 = temp1->next;
      }
      if (dir > (length / 2)) {
	temp_int = temp->start_x; temp->start_x = temp->end_x;
	temp->end_x = temp_int;
	temp_int = temp->start_y; temp->start_y = temp->end_y;
	temp->end_y = temp_int;
	temp->radius = (temp->radius < 9999.0)?(-temp->radius):999999.0;
	reverse_list( &(temp->start));
      }
      
      temp1 = temp->start; for (i=0; i< (length/2); i++) temp1 = temp1->next;
      grad_dir = temp1->grad_dir;
      sx = sx / (float)length;
      sy = sy / (float)length;
      slope = slope / (float)length;
      mag_plus = mag_plus / (float)length;
      width_plus = width_plus / (float)length;
      mag_minus = mag_minus / (float)length;
      width_minus = width_minus / (float)length;
      curvature = curvature / (float)length;
      if (temp->radius < 9999.0) {
	slope_0 = atan2 (- (float)(temp->start_x - temp->xc),
			 (float) (temp->start_y - temp->yc));
	slope_n = atan2 (- (float)(temp->end_x - temp->xc),
			 (float) (temp->end_y - temp->yc));
      }
      else {
	slope_0 = atan2 ((float) temp->start_y - temp->end_y, 
			 (float) temp->start_x - temp->end_x);
	slope_n = slope_0;
	/*fprintf(stderr, "(%d %d %f %f)", temp->start_y - temp->end_y, 
	  temp->start_x - temp->end_x, slope_0, slope_n);*/

      }
      
      /* NOTE THAT WE SAVE THE EDGES WITH THE X-Y COORDS REVERSED.
	 THUS, EACH POINT IS SPECIFIED BY (ROW, COLS) WHEREAS THE
	 WORKING XY CO-ORDS WERE ORIENTED ALONG ROW-COLS.*/

      fprintf(outfile,"(%d ",id++);
      fprintf(outfile, 
	 "(%8.3f %8.3f %8.3f %8.3f (%d %d) (%d %d) %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f (%8.3f %8.3f) %8.3f %8.3f)",
	   temp->yc,temp->xc,(temp->radius < 9999.0)?(-temp->radius):999999.0,
	      temp->fit_error/length, 
	      temp->start_y, temp->start_x, temp->end_y,temp->end_x, 
	      slope, 0.0, mag_plus,  width_plus, mag_minus, width_minus, 
	      grad_dir, curvature, sy, sx, slope_0, slope_n);
      /*if (temp->radius < 9999.0) 
	fprintf(stderr, "\n %8.3f %8.3f %8.3f", 
	temp->yc,temp->xc,(fabs(temp->radius)));*/
      
      fprintf(outfile, " (");
      temp1 = temp->start;
      while (temp1) {
	fprintf(outfile, "(%d %d) ", temp1->y, temp1->x);
	temp1 = temp1->next;
      }
      fprintf(outfile, ")) \n");
    }
    temp = temp->next;
  }
  fprintf(outfile,")");
  fclose(outfile);
}
