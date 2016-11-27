#include "header.h"

int    img_width, img_height ;     /* width and height of original image */
int    gksize, seg_focus_p;         /* size the image is expanded by */
int    exp_img_size, img_size ;    /* original and expanded image size */
int    expnd_width, expnd_height ; /* expanded width and height */
UCHAR    *image_ptr, *edge_ptr;      /*pointers to the image and edge data*/
float  b1, b2, b3, a0, a1, a2;     /* filter coefficients */
float  a01, a02, beta, alpha, tau_x, tau_y;
float    Gamma;
float  pi;
char  imgname[100], convlfilename[100];
float *convl_ptr, *dir_ptr, *slope_ptr;
int    mark_p, chain_p;
float norm_const;
float scale;



int MIN_SEG_LENGTH;
float SEG_INDEX = 200.0;
float MIN_DISTANCE = 0.5;

int  SEGMENT_P, OPERATOR = 2, MRC = 6;

float ACCPT_ERROR =  900.0;

int THRESHOLD_P, THRESH_VAL, THRESH_LENGTH;
float LOW_HYSTERESIS, HIGH_HYSTERESIS, LOW_TH, HIGH_TH;
int ATTS_FILE_P;
static int contour_id, parameter_id;


main(int argc, char **argv)
{

  struct chain *edge_chain;
  int i;

  if (argc < 3)
    {
      fprintf(stderr,  "Usage: %s img_file scale -edge 1/2 -mrc canny_mrc -lh Low Hysteresis -hh High Hysteresis -lt Length Threshold -seg segment_p -o atts_file_p\n", *argv);
      fprintf(stderr,  "       - the thresholding, the segmenting, and the edge file flags are optional.\n");
      fprintf(stderr,  "       - slope_thresh refers to the threshold on the edge profile slope\n");
      fprintf(stderr,  "       - atts_file if you want the .atts file\n");
      exit(1);
    }

  /* default parameters */
  MIN_SEG_LENGTH = 10;
  THRESHOLD_P = THRESH_VAL = THRESH_LENGTH = 0;
  LOW_HYSTERESIS = 0.0; HIGH_HYSTERESIS = 0.0;
  SEGMENT_P = 0; ATTS_FILE_P = 0;
  strcpy(convlfilename, imgname);
  
  strncpy(imgname, argv[1], 100) ;
  scale = atof(argv[2]);
  for (i=3; i < argc; i++) {
    if (strcmp(argv[i], "-edge") == 0) {OPERATOR = atoi(argv[++i]);}
    if (strcmp(argv[i], "-lh") == 0) {
      THRESHOLD_P = 1; LOW_HYSTERESIS = atof(argv[++i]);
    }
    if (strcmp(argv[i], "-hh") == 0) {
      THRESHOLD_P = 1; HIGH_HYSTERESIS = atof(argv[++i]);
    }
    /*if (strcmp(argv[i], "-st") == 0) {
      THRESHOLD_P = 1; THRESH_VAL = atoi(argv[++i]);
      }*/
    if (strcmp(argv[i], "-lt") == 0) {
      THRESHOLD_P = 1; THRESH_LENGTH = atoi(argv[++i]);
      }
    if (strcmp(argv[i], "-seg") == 0) {SEGMENT_P = atoi(argv[++i]);}
    if (strcmp(argv[i], "-o") == 0) {
      ATTS_FILE_P = 1; strcpy(convlfilename, argv[++i]);
    }
    if (strcmp(argv[i], "-mrc") == 0) {MRC = atoi(argv[++i]);}
  }
  

  read_image (imgname, &image_ptr);
  edge_chain = NULL;
  
  Gamma = scale;
  tau_x = 1.0;
  tau_y = 1.0;
  beta = 1.0 / Gamma;
  pi = 22.0/7.0;
  alpha = 0.312;


  zero_cross_filter();
  edge_chain = chain_coder();
  // if (SEGMENT_P == 1) edge_chain = segment(edge_chain);
  if (THRESHOLD_P == 1) {
    threshold_edges_length (edge_chain, (float) THRESH_LENGTH);
    threshold_edges (edge_chain, HIGH_TH);
  }
  
  
  write_chain(convlfilename, &(edge_chain));
}
