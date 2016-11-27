#ifndef _FMMISOPROPAGATION_2D_H_
#define _FMMISOPROPAGATION_2D_H_

// some global variables
extern int n;			// width
extern int p;			// height
extern double* D;
extern double* S;
extern double* W;
extern double* Q;
extern double* L;
extern double* start_points;
extern double* end_points;
extern double* H;
extern double* values;
extern int nb_iter_max;
extern int nb_start_points;
extern int nb_end_points;

typedef bool (*T_callback_insert_node)(int i, int j, int ii, int jj);

// main function
void fmmisopropagation_2d(T_callback_insert_node callback_insert_node);

#endif // _FMMISOPROPAGATION_2D_H_