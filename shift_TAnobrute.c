/* First block is taken from TA_common.c, 2nd from TA_shrink_delta, third from TA_shrink_bardelta.
All 3 were done by MWahlström. 
Any changes to the code should be accompanied by a "modif" comment.


Rest of the code is the same algorithm as shift_v2.c, where DEM algo has been replaced by TA.
Appropriate runtime limits should be looked for. 100 000 iterations were chosen for the TA like in the GWW paper.

*/



/* TO DO:
->check best_bord is updated correctly: is the best point we had previously kept? If not, need extra storage-> Looks like there are already the variables? Does the old code update them normally?
-> Store box coords appropriately-> DONE?
-> Check if box is under or over-filled-> delta seems to be doing open boxes and bardelta-> run both take worst and define fill with this.-> DONE?
-> Check where runtime cap is included?
*/




#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <string.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

typedef struct{
  double *point;
  int index;
} point; 

typedef struct{
  int key;
  double value;
} key_value;

int n_dimensions, n_points;
double *n_coords;
double **coord;
int **point_index;

int comparedim;
bool curr_fill,fill;
double *best_bord, *curr_bord;
point *optiset;
int *is_in;
double epsilon;
double globallower;
double glob_bound;
int debugg;

clock_t end, start;		/* timing variables		*/
double cput;

//Split this file into main and shift parts

int cmpdbl(const void *a, const void *b)
{
  key_value *x=(key_value*)a;
  key_value *y=(key_value*)b;
  if (x->value < y->value)
    return -1;
  else if (x->value == y->value)
    return 0;
  return 1;
}

int cmpkeyk(const void *pt1, const void *pt2)
{
  double a=(*(double **)pt1)[comparedim], b=(*(double **)pt2)[comparedim];
  if (a<b)
    return -1;
  else if (a>b)
    return 1;
  return 0;
}

void copy_point(point* dest, point* source, int dim) {
  memcpy(dest->point, source->point, n_dimensions*sizeof(double));
  dest->index = source->index;
}

void usage()
{
  fprintf(stderr, "Usage: dem_discr [dim npoints] [file]\n\nIf file not present, read from stdin. If dim, npoints not present, \nassume header '%%dim %%npoints reals' (e.g. '2 100 reals') in file.\n");
}
/****************************************
FIrst part
****************************************/

// not recommended: PRINT_ALL_UPDATES
//#define PRINT_ALL_UPDATES
// even more ridiculous!
//#define PRINT_UPDATE_CANDIDATES
//#define DISPLAY_CANDIDATES
//#define PRINT_RANGE_DATA

// for stupid speedup reasons (alternative: make it static, and take care of initialization somehow)
int *coordinate;

// we want to use C library qsort to sort.  
// I made a replacement stump
int dbl_compare(const void *a, const void *b)
{
  const double *xp=(double*)a, *yp=(double*)b;
  const double x=*xp, y=*yp;
  if (x<y)
    return -1;
  if (x>y)
    return 1;
  return 0;
} // oops: "x-y" gets cast to integer, rounded wrong

void quicksort(int left, int right, double *arr) 
{
  qsort(&arr[left], right-left+1, sizeof(double), dbl_compare);
}

double get_coord(int d, int i) 
{
  return coord[d][i];
}

int get_index_up(int d, double key)
{
  int i=(int)key*(n_coords[d]-1);
  while (coord[d][i] < key)
    i++;
  while ((i>0) && (coord[d][i-1] >= key))
    i--;
  // termination: at first coordinate which is >= key
  // (if i=0 then key==0)
  return i;
}

int get_index_down(int d, double key)
{
  int bound=n_coords[d]-1;
  int i=key*bound;
  while (coord[d][i] > key)
    i--;
  while ((i < bound) && (coord[d][i+1] <= key))
    i++;
  return i;
}

// The following two functions use an "output variable" provided by the caller
// (anything else would just mess up the C memory management)
void round_point_up(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    output[i] = get_index_up(i, point[i]);
}

void round_point_down(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    output[i] = get_index_down(i, point[i]);
}

void round_point_extradown(double *point, int *output)
{
  int i;
  for (i=0; i<n_dimensions; i++) {
    output[i] = get_index_down(i, point[i]);
    if(output[i]==0)
      output[i]=n_coords[i]-1;
  }
}

void process_coord_data(point *points, int n, int d)
{
  int i,j;
  double tmp_coords[n+2];
  int n_uniq, idx;
  double prev_coord;
  //initialise n_coords[], coord[][]
  n_dimensions=d;
  n_points=n;
  coordinate=(int*)malloc(d*sizeof(int));
  for (i=0; i<d; i++)
    coordinate[i]=i;
  n_coords = (double*)malloc(n_dimensions*sizeof(double));
  coord = (double**)malloc(n_dimensions*sizeof(double *));
  for (i=0; i<n_dimensions; i++) {
    for (j=0; j<n_points; j++)
      tmp_coords[j+1]=points[j].point[i];
    tmp_coords[0]=0.0;
    tmp_coords[n+1]=1.0;
    quicksort(1, n, tmp_coords);
    // 1. count
    n_uniq=1;
    prev_coord=0;
    for (j=1; j<=n+1; j++) { // inclusive bound: contains 1.0
      if (prev_coord==tmp_coords[j])
	continue;
      n_uniq++;
      prev_coord=tmp_coords[j];
    }
    // 2. transfer
    coord[i]=(double*)malloc(n_uniq*sizeof(double));
    idx=1;
    prev_coord=tmp_coords[0];
    coord[i][0]=prev_coord;
    for (j=1; j<=n+1; j++) {
      if (prev_coord==tmp_coords[j])
	continue;
      prev_coord=tmp_coords[j];
      coord[i][idx++] = prev_coord;
    }
    n_coords[i]=n_uniq;
//    fprintf(stderr, "Coordinates %d, %d: ", i, n_uniq);
//    for (j=0; j<n_uniq; j++)
//      fprintf(stderr, "%g ", coord[i][j]);
//    fprintf(stderr,"\n");
  }
  // finished setup for: n_coords[], coord[][]

  // next: transfer point set to into point_index
  point_index=(int**)malloc(n_points*sizeof(int *));
  for (i=0; i<n_points; i++)
    point_index[i]=(int*)malloc(n_dimensions*sizeof(int));
  for (i=0; i<n_points; i++)
    for (j=0; j<n_dimensions; j++) {
      idx=get_index_up(j, points[i].point[j]);
      if (coord[j][idx] != points[i].point[j]) {
	fprintf(stderr, "ERROR: located incorrect coordinate (%g at %d, wanted %g).\n",
		coord[j][idx], idx, points[i].point[j]);
	abort();
      }
      point_index[i][j] = idx;
    }
  // setup finished.
}

double volume(int *corner)
{
  double vol=1.0;
  int i;
  for (i=0; i<n_dimensions; i++)
    vol *= get_coord(i, corner[i]);
  return vol;
}


//is x in [0,z[ ?
int open(int *x, int *z)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    if (x[i] >= z[i])
      return 0;
  return 1;
}

//is x in [0,z] ?
int closed(int *x, int *z)
{
  int i;
  for (i=0; i<n_dimensions; i++)
    if (x[i] > z[i])
      return 0;
  return 1;
}

int count_open(int *corner)
{
  int n=n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += open(point_index[i], corner);
  return res;
}

int count_closed(int *corner)
{
  int n=n_points;
  int i;
  int res=0;
  for (i=0; i<n; i++)
    res += closed(point_index[i], corner);
  return res;
}

// group of functions for grow/"snap up"
int point_critical_dimension(int *corner, int *point)
{
  int i;
  int crit=-1;
  //  fprintf(stderr, "Point");
  //  for (i=1; i<=d; i++)
  //    fprintf(stderr, " %g", point[i]);
  for (i=0; i<n_dimensions; i++)
    if (point[i] >= corner[i]) {
      if (crit>=0) {
	//	fprintf(stderr, " double out (%d,%d)\n", crit, i);
	return -1;
      }
      else
	crit=i;
    }
  //  fprintf(stderr, " crit %d\n", crit);
  return crit;
}

// "really old" O(nd^2) time alg.
void old_grow_box(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap;
  int curr_d;
  int max_idx;
  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  /* 
  fprintf(stderr, "Growing base ");
  for (i=0; i<d; i++)
    fprintf(stderr, "%d ", corner[i]);
  fprintf(stderr, "in order ");
  for (i=0; i<d; i++)
    fprintf(stderr, "%d ", order[i]);
  fprintf(stderr, "\n");
   */

  for (i=0; i<d; i++) {
    curr_d=order[i];
    max_idx=n_coords[curr_d]-1;    //    max_val=1.0;
    /* 
    fprintf(stderr, "Direction %d (base", curr_d);
    for (j=0; j<d; j++)
      fprintf(stderr ," %d", corner[j]);
    fprintf(stderr, "), start %d\n", max_idx);
     */
    for (j=0; j<n; j++) {
      if (point_critical_dimension(corner, point_index[j])==curr_d)
	if (point_index[j][curr_d] < max_idx) {
	  max_idx=point_index[j][curr_d];
	  /* 
	  fprintf(stderr, "Because of point ");
	  for (k=0; k<d; k++)
	    fprintf(stderr,"%d ", point_index[j][k]);
	  fprintf(stderr, "bounded at %d\n", max_idx);
	   */
	}
    }
    corner[curr_d] = max_idx;
  }
}

//void grow_box_newer(int *corner)
void grow_box_randomly(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap, memo;
  int curr_d;
  int new_box[d];

  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  for (i=0; i<d; i++)
    new_box[i] = n_coords[i]-1;

  for (i=0; i<n; i++) {
    memo=-1;
    for (j=0; j<d; j++) {
      curr_d=order[j];
      k=point_index[i][curr_d];
      if (k >= corner[curr_d]) {
	if (k < new_box[curr_d]) {
	  if (memo<0)
	    memo=curr_d;
	}
	else {
	  memo=-1;
	  break;
	}
      }
    }
    if (memo >= 0)
      new_box[memo] = point_index[i][memo];
  }

#ifdef DEBUG
  if (count_open(corner) != count_open(new_box)) {
    fprintf(stderr, "ERROR: Went from %d to %d points.\n", count_open(corner), count_open(new_box));
    fprintf(stderr, "Old box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", corner[i]);
    fprintf(stderr, "\nNew box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", new_box[i]);
    fprintf(stderr, "\n");
    abort();
  }
#endif

  for (i=0; i<d; i++)
    corner[i]=new_box[i];
}

//Older version does not perform "memo" check
void grow_box_older(int *corner)
{
  int order[n_dimensions];
  int n=n_points, d=n_dimensions;
  int i,j,k, swap;
  int curr_d;
  int new_box[d];

  for (i=0; i<d; i++)
    order[i]=i;
  for (i=0; i<d; i++) {
    j = i + random()%(d-i);
    if (i != j) {
      swap=order[i];
      order[i]=order[j];
      order[j]=swap;
    }
  }
  for (i=0; i<d; i++)
    new_box[i] = n_coords[i]-1;

  for (i=0; i<n; i++) {
    for (j=0; j<d; j++) {
      curr_d=order[j];
      k=point_index[i][curr_d];
      if (k >= corner[curr_d]) {
	if (k < new_box[curr_d])
	  new_box[curr_d]=k;
	break;
      }
    }
  }
#ifdef DEBUG
  if (count_open(corner) != count_open(new_box)) {
    fprintf(stderr, "ERROR: Went from %d to %d points.\n", count_open(corner), count_open(new_box));
    fprintf(stderr, "Old box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", corner[i]);
    fprintf(stderr, "\nNew box: ");
    for (i=0; i<d; i++)
      fprintf(stderr, "%d ", new_box[i]);
    fprintf(stderr, "\n");
    abort();
  }
#endif

  for (i=0; i<d; i++)
    corner[i]=new_box[i];
}

//#define grow_box_randomly grow_box_older

// Rounds box down to borders given by its contained points.
void snap_box(int *corner)
{
  int d=n_dimensions;
  int n=n_points;
  int max_idx[d];
  int i,j;
  for (i=0; i<d; i++)
    max_idx[i]=-1;
  for (i=0; i<n; i++)
    if (closed(point_index[i], corner))
      for (j=0; j<d; j++)
	if (point_index[i][j] > max_idx[j])
	  max_idx[j]=point_index[i][j];
  for (i=0; i<d; i++)
    corner[i] = (max_idx[i] < 0) ? 0 : max_idx[i];
}



//calculates delta(x)
double get_delta(int *x)
{ 
  int i, op;
  double vol, delta;
  int n=n_points;
  
  //calculates no of points in [0,x[
  op = count_open(x);

  //calcualtes volume of box generated by x
  vol = volume(x);
  
  //calculates delta
  delta = vol - (double)op/n;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, op, (double)op/n, delta);  
  return delta;
}

//calculates bar(delta)(x)
double get_bar_delta(int *x)
{ 
  int i, cl;
  int n=n_points;
  double vol, bdelta;
  
  //calculates no of points in [0,x]
  cl = count_closed(x);
  
  //calcualtes volume of box generated by x
  vol = volume(x);
  
  //calculates bar(delta)
  bdelta = (double)cl/n - vol;
  //  fprintf(stderr, "Stand.: vol %g pts %d -> %g error %g\n",
  //	  vol, cl, (double)cl/n, bdelta);  
  return bdelta;
}

//Generate random search point xc
//Fills coordinates (indexes, not numbers) into its three arguments
void generate_xc(int *xn_plus, int *xn_minus, int *xn_extraminus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_up(xn, xn_plus);
  round_point_down(xn, xn_minus);
  round_point_extradown(xn, xn_extraminus);		
}

void generate_xc_delta(int *xn_plus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_up(xn, xn_plus);
}

void generate_xc_bardelta(int *xn_minus, int *xn_extraminus)
{
  int j, d=n_dimensions;
  double xn[d];
  double temp;
  for(j=0; j<d; j++)
    {
      temp=(double)((double)rand()/RAND_MAX);
      xn[j]=pow(temp,(double)((double)1/(double)d));
    }
  round_point_down(xn, xn_minus);
  round_point_extradown(xn, xn_extraminus);		
}


//Generate a random neighbor of xc
//k[i] == range radius in component i, mc = number of components to change
//the three xn_* variables are filled with indexes
void generate_neighbor (int *xn_plus_index, int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(xn, xn_plus_index);
  round_point_down(xn, xn_minus_index);
  round_point_extradown(xn, xn_extraminus_index);
}

void generate_neighbor_delta(int *xn_plus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_up(xn, xn_plus_index);
}

void generate_neighbor_bardelta(int *xn_minus_index, int *xn_extraminus_index, 
			int *xc_index, int *k, int mc)
{
  int i, j, q, d=n_dimensions;
  double temp, upper_bound, lower_bound;
  double xn[d];

  //First copy the values of the current search point 
  for(j=0; j<d; j++)
    {
      xn[j] = coord[j][xc_index[j]];
    }

  // find mc different coordinates to be changed
  for (j=0; j<mc; j++) {
    i = j + random()%(d-j);
    if (i != j) {
      q = coordinate[j];
      coordinate[j]=coordinate[i];
      coordinate[i] = q;
    }
  }
  
  //set lower and upper bound to the box from which the random neighbor will be sampled
  for(j=0; j<mc; j++){ 	
    if (xc_index[coordinate[j]]-k[coordinate[j]]>=0) 
      lower_bound = coord[coordinate[j]][xc_index[coordinate[j]]-k[coordinate[j]]];
    else 
      lower_bound=0.0;
    if (xc_index[coordinate[j]]+k[coordinate[j]] < n_coords[coordinate[j]])
      upper_bound = coord[coordinate[j]][xc_index[coordinate[j]]+k[coordinate[j]]];
    else 
      upper_bound=1.0;
    
    //draw a random number in [0,1]
    temp=(double)((double)rand()/RAND_MAX);
    temp=(pow(upper_bound,d)-pow(lower_bound,d))*temp + pow(lower_bound,d);
    xn[coordinate[j]]=pow(temp,(double)((double)1/(double)d));
  }

  round_point_down(xn, xn_minus_index);
  round_point_extradown(xn, xn_extraminus_index);
}



/******************************************************
Second part
******************************************************/

#ifndef MC
#define MC 2
#endif

// not recommended: PRINT_ALL_UPDATES
//#define PRINT_ALL_UPDATES
// even more ridiculous!
//#define PRINT_UPDATE_CANDIDATES
//#define DISPLAY_CANDIDATES
//#define PRINT_RANGE_DATA

// use THRESH_REPEAT=1 in "production code" (it only helps to denoise testing)
#define THRESH_REPEAT 1

int k_div=0;// "0" means default "4 or 8" setup, other value (from main(), e.g. -k 16) overrides

#define I_TILDE 316    // thresholds to be calculated (sqrt(iterations)), default value 100k
int mc=MC;             //nbr of coordinates to be changed, default value
int i_tilde=I_TILDE;
#define TRIALS 10      //nbr of runs (mean and max will be calculated), default value
int trials=TRIALS;

#ifndef VERSION
#define VERSION 0
#endif

// global variables to store info about pointset


// global variables to store info about worst box (also "private")
double real_max_discr=0;
int real_when=0, when=0;
int current_iteration;


//Computes the best of the rounded points -- basic version
//Constant neighbourhood size and mc-values.
//Does not split the search.
//Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_delta(int *xn_plus)
{
  double fxc;  
  int j, d=n_dimensions;
#if (VERSION == 1) || (VERSION == 2)
  int xn_plus_grow[d];
#endif

  // Growing, shrinking.
#if (VERSION == 3)
  // SNAPUPDATE
  grow_box_randomly(xn_plus);
#elif (VERSION == 1) || (VERSION == 2)
  // Grower, shrinker that copy the point
  for (j=0; j<d; j++)
    xn_plus_grow[j]=xn_plus[j];
  grow_box_randomly(xn_plus_grow);
#endif

  // For version 1: "private update"
#if (VERSION == 1)
  fxc = get_delta(xn_plus_grow);
#ifdef PRINT_UPDATE_CANDIDATES
  fprintf(stderr, "PRIVATE candidate %g (vs %g)\n",
	  fxc, real_max_discr);
#endif
  if (fxc > real_max_discr) {
    real_max_discr = fxc;
    real_when = current_iteration;
#ifdef PRINT_ALL_UPDATES
    fprintf(stderr, "Secret update at %d to %g\n",
	    current_iteration, fxc);
#endif
  }
#endif
  // Ends "private update" block

  // Now, create the official numbers.
#if (VERSION == 2)
  // official update from modified points
  fxc = get_delta(xn_plus_grow);
#else
  // versions 0,3 both compute from the point now given by xn_*
  // version 1 officially reports this as well
  fxc = get_delta(xn_plus);
#endif

  // In delta_only mode: don't copy

  return fxc;
}



double oldmain(point *pointset, int n, int d) 
{
  int k[d], start[d];
  
  int i, j, p, t;           // loop variables
  
  double thresh[i_tilde];    //Thresholdsequence
  double T;                            //current Threshold
  
  double fxc;
  int xc_index[d], xn_plus_index[d];
  int xn_best_index[d];    //Indices of current point, neighbour
  double xglobal[trials+1][d], xbest[d];
  double current, global[trials+1], best, mean;  //current and global best values
  
  int outerloop=i_tilde, innerloop=i_tilde;     
  
  int anzahl=0;
  int switches[trials+1]; 
  int global_switches[trials+1];    
  
  //Get pointset from external file
  FILE *datei_ptr=stderr;
  //fprintf(datei_ptr,"GLP-Menge %d %d  ",d,n);
  
  //Sort the grid points, setup global variables
  process_coord_data(pointset, n, d);
  
  //Algorithm starts here
  for(t=1;t<=trials;t++)
    { //Initialization
      //fprintf(stderr, "Trial %d/%d\n", t, trials);

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2;

      //Initialize iteration count
      current_iteration=0;

      //Generate threshold sequence
      for(i=1;i<=outerloop;i++){
	
	current_iteration++;
	//Update k-value
	  for (j=0; j<d; j++) {
	    k[j] = start[j]*(((double)outerloop-current_iteration)/(outerloop)) +
	      1*((double)current_iteration/(outerloop));
		  //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
	  }

        //Update mc-value
	  mc=2+(int)(current_iteration/outerloop*(d-2));


	//generation of random point xc
	generate_xc_delta(xc_index); 
	
	//(Possibly) Snaps the point upwards and computes the fitness
	current = best_of_rounded_delta(xc_index);
    
	//draw a neighbour of xc
	generate_neighbor_delta(xn_plus_index, xc_index, k, mc);
	
	//Compute the threshold
	fxc=best_of_rounded_delta(xn_plus_index);
	thresh[i]=0.0-fabs(fxc-current);
      }	
  
      //sort the thresholds in increasing order
      quicksort(1,outerloop,thresh);
  

      switches[t]=0;
      global_switches[t]=0;
      current=0;
      global[t]=0;
      when=0;
      real_when=0;
      real_max_discr=0;

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));


      //draw a random initial point 
      generate_xc_delta(xc_index);
   
      //(Possibly) Snap and compute the best of the rounded points and update current value
      current = best_of_rounded_delta(xc_index);
      
      global[t] = current;

      current_iteration=0;
      for(i=1;i<=outerloop;i++)
	{
	  T=thresh[i];
	  
	  for(p=1;p<=innerloop;p++)
	    {
	      current_iteration++;
	      
	      //Update k-value
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, "Snapshot: range ");
#endif
	      for (j=0; j<d; j++) {
		k[j] = start[j]*(((double)innerloop*outerloop-current_iteration)/(innerloop*outerloop)) +
		  1*((double)current_iteration/(innerloop*outerloop));
		  //		k[j]=(int)(start[j]-(int)(current_iteration/(innerloop*outerloop)*(start[j]-1)));
#ifdef PRINT_RANGE_DATA
		if (p==1)
		  fprintf(stderr, "%d ", k[j]);
#endif
	      }

	      //Update mc-value
	      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, " threshold %g mc %d\n", T, mc);
#endif
	      //mc=2;

	      //Get random neighbor
	      generate_neighbor_delta(xn_plus_index, xc_index,k,mc);
#ifdef DISPLAY_CANDIDATES
	      fprintf(stderr, "Old: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xc_index[j]);	    
	      fprintf(stderr, "\nPlus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_plus_index[j]);
	      fprintf(stderr, "\n");
#endif

	      //(Possibly) Snap the points and compute the best of the rounded points 
	      fxc = best_of_rounded_delta(xn_plus_index);
#ifdef PRINT_UPDATE_CANDIDATES
	      fprintf(stderr, "Iter. %d candidate %10g (vs %10g best %10g) -- ",
		      current_iteration, fxc, current, global[t]);
#endif
	      //Global update if necessary
	      if(fxc>global[t]){
		global_switches[t]++;
		global[t]=fxc;
		when=current_iteration;
		
		//MODIF: looks like boxes weren't being updated-> code at the end was a stump
		for (j=0;j<d;j++){
			xglobal[t][j]= xn_plus_index[j]; // IS THIS CORRECT? IT'S WHAT WAS BEING PRINTED...
		}
		
		
		
		
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "global ");
#endif
#ifdef PRINT_ALL_UPDATES
		fprintf(stderr, "%g at %d :", fxc, current_iteration);
		for (j=0; j<d; j++)
		  fprintf(stderr, " %d", xn_plus_index[j]);
		fprintf(stderr, "\n");
#endif
	      }
	      //Update of current best value if necessary
	      if(fxc-current>=T){
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "update\n");
#endif
		switches[t]++;
		current=fxc;
		for(j=0; j<d; j++){
		  xc_index[j]=xn_plus_index[j];
		}
	      }
#ifdef PRINT_UPDATE_CANDIDATES
	      else {
		fprintf(stderr, "skip\n");
	      }
#endif
	    }//innerloop
	}//outerloop
      if (real_max_discr > global[t]) {
	global[t] = real_max_discr;
	when = real_when;
	//	fprintf(stderr, "Max value subsumed\n");
      }
      //fprintf(stderr, "Result %g at %d\n", global[t], when);
      //fprintf(stdout, "%g\n", global[t]); // To simplify post-execution bookkeeping
    }//trials
  
  
  //best calculated value 
  best=global[1];
  for(j=0; j<d; j++) xbest[j]=xglobal[1][j];
  for(t=2;t<=trials;t++)
    {
      if(global[t]>best)
	{ 
	  best=global[t];
	  for(j=0; j<d; j++) xbest[j]=xglobal[t][j];
	}
    }
  for (j=0;j<d;j++){ //modif here
	  best_bord[j]=get_coord(j, xbest[j]);
  }
  
  for(t=1;t<=trials;t++)
    {
      if(global[t]==best) anzahl++;
    }
  //fprintf(datei_ptr,"best %e  ",best);
  // for(j=0; j<d; j++)  fprintf(datei_ptr,"xbest %d coo  %e\n", j,xbest[j]);
  
  //delta or bar(delta) causing best value?
  //if(best==fabs(delta(xbest,GLP))) fprintf(datei_ptr,"delta\n");
  //else fprintf(datei_ptr,"bar_delta\n");
  
  //calculation of mean value
  mean=0;
  for(t=1;t<=trials;t++) mean=mean+global[t];
  mean=mean/trials;
  //fprintf(datei_ptr,"mean %e  ",mean);
  //fprintf(datei_ptr,"lower_bound %e\n",lower_bound);
  //fprintf(datei_ptr,"upper_bound %e\n",upper_bound);
  
  //  fprintf(datei_ptr,"Anzahl der Iterationen: %d  ",iteration_count);
  //fprintf(datei_ptr,"Wert von k: %d\n",k);
  // fprintf(datei_ptr,"Wert von Extraminus: %d\n",extraminus);
  //fprintf(datei_ptr,"Anzahl best: %d\n",anzahl);
  // for(i=1;i<=outerloop;i++) fprintf(datei_ptr,"Thresh %d = %e\n",i,thresh[i]);  
  
  // for(t=1;t<=trials;t++) { 
  //fprintf(datei_ptr,"Anzahl switches in Runde %d: %d\n",t,switches[t]);
  //fprintf(datei_ptr,"Anzahl global_switches in Runde %d: %d\n",t,global_switches[t]);
  //}
  
  //MODIF: there were no frees
  free(coordinate);
  free(n_coords);
  for (i=0;i<n_dimensions;i++)
	  free(coord[i]);
  free(coord);
  for (i=0;i<n_points;i++)
	  free(point_index[i]);
  free(point_index);
  
  return best;
  
}


/*****************************************
Third part
*****************************************/



//Computes the best of the rounded points -- basic version
//Constant neighbourhood size and mc-values.
//Does not split the search.
//Copies the appropriate "thing" into xc_index (output variable)
double best_of_rounded_bardelta(int *xn_minus, int *xn_extraminus, int *xc_index)
{
   double fxn_minus;
  double fxn_extraminus;
  double fxc;  
  int j, d=n_dimensions;
  int use_extraminus=0;
#if (VERSION == 1) || (VERSION == 2)
  int xn_minus_snap[d], xn_extraminus_snap[d];
#endif
  for (j=0; j<d; j++)
    if (xn_minus[j] != xn_extraminus[j]) {
      use_extraminus=1;
      break;
    }

  // Growing, shrinking.
#if (VERSION == 3)
  // SNAPUPDATE
  snap_box(xn_minus);
  if (use_extraminus)
    snap_box(xn_extraminus);
#elif (VERSION == 1) || (VERSION == 2)
  // Grower, shrinker that copy the point
  for (j=0; j<d; j++)
    xn_minus_snap[j]=xn_minus[j];
  snap_box(xn_minus_snap);
  if (use_extraminus) {
    for (j=0; j<d; j++)
      xn_extraminus_snap[j]=xn_extraminus[j];
    snap_box(xn_extraminus_snap);
  }
#endif

  // For version 1: "private update"
#if (VERSION == 1)
  fxc = get_bar_delta(xn_minus_snap);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus_snap);
    fxc=max(fxc, fxn_extraminus);
  }
#ifdef PRINT_UPDATE_CANDIDATES
  fprintf(stderr, "PRIVATE candidate %g (vs %g)\n",
	  fxc, real_max_discr);
#endif
  if (fxc > real_max_discr) {
    real_max_discr = fxc;
    real_when = current_iteration;
#ifdef PRINT_ALL_UPDATES
    fprintf(stderr, "Secret update at %d to %g\n",
	    current_iteration, fxc);
#endif
  }
#endif
  // Ends "private update" block

  // Now, create the official numbers.
#if (VERSION == 2)
  // official update from modified points
  fxc = get_bar_delta(xn_minus_snap);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus_snap);
    fxc=max(fxc, fxn_extraminus);
  }
#else
  // versions 0,3 both compute from the point now given by xn_*
  // version 1 officially reports this as well
  fxc = get_bar_delta(xn_minus);
  if (use_extraminus) {
    fxn_extraminus=get_bar_delta(xn_extraminus);
    fxc=max(fxc, fxn_extraminus);
  }
#endif

  // Remains only to copy the winning point to output variable xc_index.
  if (use_extraminus && (fxn_extraminus >= fxc)) {
    for (j=0; j<d; j++)
      xc_index[j] = xn_extraminus[j];
  }
  else { 
    for(j=0; j<d; j++) 
      xc_index[j]= xn_minus[j];
  }

  return fxc;
}



double oldmainbar(point *pointset, int n, int d)
{
  int k[d], start[d];
  
  int i, j, p, t;           // loop variables
  
  double thresh[i_tilde];    //Thresholdsequence
  double T;                            //current Threshold
  
  double fxc;
  int xc_index[d], xn_minus_index[d], xn_extraminus_index[d];
  int xn_best_index[d];    //Indices of current point, neighbour
  double xglobal[trials+1][d], xbest[d];
  double current, global[trials+1], best, mean;  //current and global best values
  
  int outerloop=i_tilde, innerloop=i_tilde;     
  
  int anzahl=0;
  int switches[trials+1]; 
  int global_switches[trials+1];    
  
  //Get pointset from external file
  FILE *datei_ptr=stderr;
  //fprintf(datei_ptr,"GLP-Menge %d %d  ",d,n);
  
  //Sort the grid points, setup global variables
  process_coord_data(pointset, n, d);
  
  //Algorithm starts here
  for(t=1;t<=trials;t++)
    { //Initialization
      //fprintf(stderr, "Trial %d/%d\n", t, trials);

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2;

      //Initialize iteration count
      current_iteration=0;

      //Generate threshold sequence   (only once)
      //      fprintf(stderr, "Generating threshold\n");
      for(i=1;i<=outerloop;i++){
	
	current_iteration++;
	//Update k-value
	  for (j=0; j<d; j++) {
	    k[j] = start[j]*(((double)outerloop-current_iteration)/(outerloop)) +
	      1*((double)current_iteration/(outerloop));
		  //	    k[j]=start[j] - (int)((3.0/4)*(current_iteration/outerloop)*(start[j]-1));
	  }

        //Update mc-value
	  mc=2+(int)(current_iteration/outerloop*(d-2));


	//generation of random point xc
	generate_xc_bardelta(xn_minus_index, xn_extraminus_index); 
	
	//(Possibly) Snap the points and compute the largest of the rounded values 
	current = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);
    
	//draw a neighbour of xc
	generate_neighbor_bardelta(xn_minus_index, xn_extraminus_index, xc_index, k, mc);
	
	//Compute the threshold
	fxc=best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);
	thresh[i]=0.0-fabs(fxc-current);
      }	
  
      //sort the thresholds in increasing order
      quicksort(1,outerloop,thresh);
  

      switches[t]=0;
      global_switches[t]=0;
      current=0;
      global[t]=0;
      when=0;
      real_when=0;
      real_max_discr=0;

      //Initialize k-value
      for (j=0; j<d; j++) {
	start[j]=(int)((n_coords[j]-1)/2);
      }
      //Initialize mc-value
      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));


      //draw a random initial point 
      generate_xc_bardelta(xn_minus_index, xn_extraminus_index);       
   
      //(Possibly) Snap and compute the best of the rounded points and update current value
      current = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xc_index);
      
      global[t] = current;

      current_iteration=0;
      for(i=1;i<=outerloop;i++)
	{
	  T=thresh[i];
	  
	  for(p=1;p<=innerloop;p++)
	    {
	      current_iteration++;
	      
	      //Update k-value
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, "Snapshot: range ");
#endif
	      for (j=0; j<d; j++) {
		k[j] = start[j]*(((double)innerloop*outerloop-current_iteration)/(innerloop*outerloop)) +
		  1*((double)current_iteration/(innerloop*outerloop));
		  //		k[j]=(int)(start[j]-(int)(current_iteration/(innerloop*outerloop)*(start[j]-1)));
#ifdef PRINT_RANGE_DATA
		if (p==1)
		  fprintf(stderr, "%d ", k[j]);
#endif
	      }

	      //Update mc-value
	      mc=2+(int)(current_iteration/(innerloop*outerloop)*(d-2));
#ifdef PRINT_RANGE_DATA
	      if (p==1)
		fprintf(stderr, " threshold %g mc %d\n", T, mc);
#endif
	      //mc=2;

	      //Get random neighbor
	      generate_neighbor_bardelta(xn_minus_index, xn_extraminus_index, xc_index,k,mc);
#ifdef DISPLAY_CANDIDATES
	      fprintf(stderr, "Old: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xc_index[j]);	    
	      fprintf(stderr, "\nMinus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_minus_index[j]);
	      fprintf(stderr, "\nXMinus: ");
	      for (j=0; j<d; j++)
		fprintf(stderr, "%d ", xn_extraminus_index[j]);
	      fprintf(stderr, "\n");
#endif

	      //(Possibly) Snap the points and compute the best of the rounded points 
	      fxc = best_of_rounded_bardelta(xn_minus_index, xn_extraminus_index, xn_best_index);
#ifdef PRINT_UPDATE_CANDIDATES
	      fprintf(stderr, "Iter. %d candidate %10g (vs %10g best %10g) -- ",
		      current_iteration, fxc, current, global[t]);
#endif
	      //Global update if necessary
	      if(fxc>global[t]){
		global_switches[t]++;
		global[t]=fxc;
		when=current_iteration;
		
		//MODIF: looks like boxes weren't being updated-> code at the end was a stump
		for (j=0;j<d;j++){
			xglobal[t][j]= xn_best_index[j]; // IS THIS CORRECT? IT'S WHAT WAS BEING PRINTED...
		}
		
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "global ");
#endif
#ifdef PRINT_ALL_UPDATES
		fprintf(stderr, "%g at %d :", fxc, current_iteration);
		for (j=0; j<d; j++)
		  fprintf(stderr, " %d", xn_best_index[j]);
		fprintf(stderr, "\n");
#endif
	      }
	      //Update of current best value if necessary
	      if(fxc-current>=T){
#ifdef PRINT_UPDATE_CANDIDATES
		fprintf(stderr, "update\n");
#endif
		switches[t]++;
		current=fxc;
		for(j=0; j<d; j++){
		  xc_index[j]=xn_best_index[j];
		}
	      }
#ifdef PRINT_UPDATE_CANDIDATES
	      else {
		fprintf(stderr, "skip\n");
	      }
#endif
	    }//innerloop
	}//outerloop
      if (real_max_discr > global[t]) {
	global[t] = real_max_discr;
	when = real_when;
	//	fprintf(stderr, "Max value subsumed\n");
      }
      //fprintf(stderr, "Result %g at %d\n", global[t], when);
      //fprintf(stdout, "%g\n", global[t]); // To simplify post-execution bookkeeping
    }//trials
  
  
  //best calculated value 
  best=global[1];
  for(j=0; j<d; j++) xbest[j]=xglobal[1][j];
  for(t=2;t<=trials;t++)
    {
      if(global[t]>best)
	{ 
	  best=global[t];
	  for(j=0; j<d; j++) xbest[j]=xglobal[t][j];
	}
    }
  
  for (j=0;j<d;j++){ //modif here
	  best_bord[j]=get_coord(j, xbest[j]);
  }
  
  for(t=1;t<=trials;t++)
    {
      if(global[t]==best) anzahl++;
    }
  //fprintf(datei_ptr,"best %e  ",best);
  // for(j=0; j<d; j++)  fprintf(datei_ptr,"xbest %d coo  %e\n", j,xbest[j]);
  
  //delta or bar(delta) causing best value?
  //if(best==fabs(delta(xbest,GLP))) fprintf(datei_ptr,"delta\n");
  //else fprintf(datei_ptr,"bar_delta\n");
  
  //calculation of mean value
  mean=0;
  for(t=1;t<=trials;t++) mean=mean+global[t];
  mean=mean/trials;
  //fprintf(datei_ptr,"mean %e  ",mean);
  //fprintf(datei_ptr,"lower_bound %e\n",lower_bound);
  //fprintf(datei_ptr,"upper_bound %e\n",upper_bound);
  
  //  fprintf(datei_ptr,"Anzahl der Iterationen: %d  ",iteration_count);
  //fprintf(datei_ptr,"Wert von k: %d\n",k);
  // fprintf(datei_ptr,"Wert von Extraminus: %d\n",extraminus);
  //fprintf(datei_ptr,"Anzahl best: %d\n",anzahl);
  // for(i=1;i<=outerloop;i++) fprintf(datei_ptr,"Thresh %d = %e\n",i,thresh[i]);  
  
  // for(t=1;t<=trials;t++) { 
  //fprintf(datei_ptr,"Anzahl switches in Runde %d: %d\n",t,switches[t]);
  //fprintf(datei_ptr,"Anzahl global_switches in Runde %d: %d\n",t,global_switches[t]);
  //}
  
  
  //MODIF: there were no frees at any point???
  free(coordinate);
  free(n_coords);
  for (i=0;i<n_dimensions;i++)
	  free(coord[i]);
  free(coord);
  for (i=0;i<n_points;i++)
	  free(point_index[i]);
  free(point_index);
  
  return best;
  
}



/********************************
ACTUAL SHIFT HEURISTIC
********************************/

int find_point(key_value **copysorted,double temp_coord, int cho, int npoints){
	int a,b,mid;
	a=0;
	b=npoints-1;
	bool found;
	found=false;
	mid=(a+b)/2;
	while (b-a>1 && !found){
		if (copysorted[cho][mid].value == temp_coord){
			return(mid);
		}
		else if (copysorted[cho][mid].value < temp_coord){
			a=mid;
			mid=(a+b)/2;
		}
		else {
			b=mid;
			mid=(a+b)/2;
		}
			 
	}
	if (b-a==1){ // This only works because we KNOW the box should be given by some coord. With double precision there might be some mistake in the volume calc-> inexact coord. We take the closest one.
		if (copysorted[cho][b].value-temp_coord> temp_coord-copysorted[cho][a].value)
			return b;
		else
			return a;
	}
	return -1;
}


// ADD THIS AFTER EACH DISCRE CALC!!
void replace(double **pointset, double **Orig_pointset, int kpoints, int npoints){
	int i,j;
	for(i=0;i<npoints;i++)
		is_in[i]=-1;
	for (i=0;i<kpoints;i++){
		for (j=0; j<npoints;j++){
			if (fabs(pointset[i][0]-Orig_pointset[j][0])<epsilon/2){
				is_in[j]=i;
			}
		}
	}
	return;
	
}



double shift(int npoints, int kpoints, int dim, point *Orig_pointset)
{
  
  // SHuffle our point array, shuffle directly in Orig_pointset. Do this in main?
  
  
  int i, j, h, a, b, c, d, f, g, search_count, temp_bordpoint, chosen_dim, nb_natural, nb_brute, nb_runs, tempo, curr_dim, actu_dim, index;
  double upper,lower, uppero,upperc;
  upper=1.0;
  int nb_calc=0;
  
	
// Sorting the points in each dim
  key_value **copysorted=(key_value**)malloc(dim*sizeof(key_value*));
  point *subset;
  point *temp_subset;
  int **orderings=(int**)malloc(dim*sizeof(int*));// The ordering[i][j]=h means that point j is h-th in the ordering in dimension i.
  int **revorderings=(int**)malloc(dim*sizeof(int*)); // revordering[i][j]=h means that the j-th point in dimension i is point h.
  for (i=0;i<dim;i++){
	  copysorted[i]=(key_value*)malloc(npoints*sizeof(key_value));
	  orderings[i]=(int*)malloc(npoints*sizeof(int));
	  revorderings[i]=(int*)malloc(npoints*sizeof(int));
	  for (j=0;j<npoints;j++) {
		  copysorted[i][j].value=Orig_pointset[j].point[i]; // Warning dimensions switched
      copysorted[i][j].key=j;
    }
	  qsort(&(copysorted[i][0]), npoints, sizeof(key_value), cmpdbl);
	  for (j=0;j<npoints;j++){// Need points in general position for this
      h = copysorted[i][j].key;
      orderings[i][h]=j;
      revorderings[i][j]=h;
	  } 
  }
  is_in =(int*)malloc(npoints*sizeof(int));
  for (i=0;i<npoints;i++){
	  if (i<kpoints)
		is_in[i]=i;
	  else
		is_in[i]=-1;
	}
	subset=(point*)malloc(kpoints*sizeof(point)); // INTRODUCE kpoints
	temp_subset=(point*)malloc(kpoints*sizeof(point));
	for (i=0;i<kpoints;i++){// Our current point set and create a future tep_subset for the possible changes.
		subset[i].point=(double*)malloc(dim*sizeof(double));
		temp_subset[i].point=(double*)malloc(dim*sizeof(double));
    copy_point(&subset[i],&Orig_pointset[i], dim);
    copy_point(&temp_subset[i],&Orig_pointset[i], dim);
	}
	fprintf(stderr, "Sorted points\n");
	//upper = oydiscr(temp_subset, dim, kpoints,&lower); // CHANGE HERE
	uppero=oldmain(temp_subset,kpoints,dim);
	upperc=oldmainbar(temp_subset,kpoints,dim);
	if (uppero>upperc){
		upper=uppero;
		fill=false;
	}
	else{
		upper=upperc;
		fill=true;
	}
	fprintf(stderr,"Init:%f",upper);
	glob_bound=0;
	
	//lower=0.0;
	nb_calc+=1;
	//globallower=0.0;
	//replace(subset,Orig_pointset,kpoints,npoints); 
	double curr_disc, curr_disco,curr_discc;
	bool chosen,boom,problem;
	nb_natural=0;nb_brute=0;
	double *top_bord;
	top_bord=(double*)malloc(dim*sizeof(double));
	//memcpy(top_bord,best_bord,dim*sizeof(double));
	for (i=0;i<dim;i++)
		top_bord[i]=best_bord[i];
	bool insidei, insidej;
	// Following initialisations shouldn't be necessary?
	temp_bordpoint=-1;
	curr_disc=1.0;
	nb_runs=1000; // Tweak it here directly (useless for the moment)
	problem=false;
	int *list_points;
	list_points=(int*)malloc(dim*sizeof(int));
	//
	for (b=0; b<nb_runs;b++){
		chosen_dim=rand() %dim;
		chosen=false;
		for (i=0;i<dim;i++){
			temp_bordpoint=find_point(copysorted,top_bord[i],i,npoints);
			c=revorderings[i][temp_bordpoint];
			if (is_in[c]< -0.5){
				problem=false;
				index=temp_bordpoint;
				while(index >=0 && is_in[revorderings[i][index]]< -0.5 )
					index--;
				if (index==-1){
					problem=true;
					break;
				}
				c=revorderings[i][index];
			}
			list_points[i]=c;
		}
		search_count=0;
		actu_dim;
		while (!chosen && search_count<npoints && !problem){ // bound on search_count could be improved
			curr_dim=0;
			if (fill){
				g=0;
				while (!chosen && curr_dim<dim){
					actu_dim=(curr_dim+chosen_dim)%dim;
					c=orderings[actu_dim][list_points[actu_dim]]; // The position of the point we want to replace. Could maybe pre-define a table to avoid recomputing this every time
					if (c+search_count>=npoints){
						g+=1;
						curr_dim+=1;
						if (g==dim) // We're stuck, no need to go further
							break;
						else
							continue;
					}
					d=revorderings[actu_dim][c+search_count]; // Candidate for replacement
					if (is_in[d]< -0.5){// The point was not already in the set
						c=list_points[actu_dim]; // Before we had the position and not the point number
						for (i=0;i<kpoints;i++)
              copy_point(&temp_subset[i],&subset[i],dim);
						tempo=is_in[c];
            copy_point(&temp_subset[tempo],&Orig_pointset[d],dim);
						//curr_disc = oydiscr(temp_subset, dim, kpoints, &lower); //CHANGE HERE
						curr_disco=oldmain(temp_subset,kpoints,dim);
						curr_discc=oldmainbar(temp_subset,kpoints,dim);
						if (curr_disco>curr_discc){
							curr_disc=curr_disco;
							curr_fill=false;
						}
						else{
							curr_disc=curr_discc;
							curr_fill=true;
						}
						glob_bound=0;
						//lower=0.0;
						//globallower=0.0;
						nb_calc+=1;
						if (curr_disc<upper) {// Our replacement is good
							chosen=true;
							nb_natural+=1;
              copy_point(&temp_subset[tempo],&Orig_pointset[d],dim);
							is_in[d]=tempo;
							is_in[c]=-1;
							fill=curr_fill;
							memcpy(top_bord,best_bord,dim*sizeof(double));
							upper=curr_disc;
							fprintf(stderr,"New:%lf",upper);
						}
						//BUG WAS HERE, WITH AN ELSE MODIFYING OUR SUBSET FOR NO REASON
						}
				curr_dim+=1;	
					}
			}
			
			else{ // Now underfilled box
				g=0;
				while (!chosen && curr_dim<dim){
					actu_dim=(curr_dim+chosen_dim)%dim;
					c=orderings[actu_dim][list_points[actu_dim]]; // The point we want to replace
					if (c-search_count<0){
						g+=1;
						curr_dim+=1;
						if (g==dim) // We're stuck
							break;
						else
							continue;
					}
					d=revorderings[actu_dim][c-search_count];
					if (is_in[d]< -0.5){// The point was not already in the set
						c=list_points[actu_dim]; // Before we had the position and not the point number
						for (i=0;i<kpoints;i++)
              copy_point(&temp_subset[tempo],&Orig_pointset[d],dim);
						tempo=is_in[c];							
            copy_point(&temp_subset[tempo],&Orig_pointset[d],dim);
						//curr_disc = oydiscr(temp_subset, dim, kpoints, &lower); // CHANGE HERE
						curr_disco=oldmain(temp_subset,kpoints,dim);
						curr_discc=oldmainbar(temp_subset,kpoints,dim);
						if (curr_disco>curr_discc){
							curr_disc=curr_disco;
							curr_fill=false;
						}
						else{
							curr_disc=curr_discc;
							curr_fill=true;
						}
						glob_bound=0;
						//lower=0.0;
						globallower=0.0;
						nb_calc+=1;
						if (curr_disc<upper) {// Our replacement is good
							chosen=true;
							nb_natural+=1;
              copy_point(&temp_subset[tempo],&Orig_pointset[d],dim);
							is_in[d]=tempo;
							is_in[c]=-1;
							
							fill=curr_fill;
							memcpy(top_bord,best_bord,dim*sizeof(double));
							
							upper=curr_disc;
							fprintf(stderr,"New2:%lf",upper);
						}
					}
						curr_dim+=1;
				}
			}
		search_count+=1; // We've gone thourgh a full dim rotation or found a point	
		}
	if (!chosen)
		break;
	}
	fprintf(stderr,"Nb calcu:%d, Final discre:%lf\n",nb_calc,upper);
	for (i=0;i<kpoints;i++)
    copy_point(&optiset[i],&subset[i],dim);
		
	for (i=0;i<dim;i++){
		free(copysorted[i]);
		free(orderings[i]);
		free(revorderings[i]);
	}
	for (i=0;i<kpoints;i++){
		free(subset[i].point);
		free(temp_subset[i].point);
	}
	free(copysorted);
	free(subset);
	free(temp_subset);
	free(orderings);
	free(revorderings);
	free(is_in);
	free(list_points);
	
	free(top_bord);
	
  fprintf(stderr,"Natural: %d, Brute: %d, Discr: %lf",nb_natural,nb_brute,upper);
  return upper;

  
}

int main(int argc, char **argv)
{
  int dim, npoints,i,j,h,kpoints;
  int chosen_dim;
  FILE *pointfile;
  double upper,lower;
  int nb_tries;
  bool fill; //Tracks if the worst box was under or over filled
  point *Orig_pointset;
  FILE *random;
  unsigned int seed;
  srand(1);
  pointfile = fopen(argv[1], "r");
  /*int test;
  test=fscanf(pointfile,"%d %d %d", &dim, &npoints, &kpoints);
  if (test!=3)
	  exit(EXIT_FAILURE);
  fprintf(stderr, "\n \n Reading dim %d npoints %d kpoints %d\n", dim, npoints,kpoints);*/
  
  dim=atoi(argv[2]);
  npoints=atoi(argv[3]);
  kpoints=atoi(argv[4]);
  fprintf(stderr, "Reading dim %d npoints %d kpoints %d ", dim, npoints,kpoints);
  
  Orig_pointset = (point*)malloc(npoints*sizeof(point));
  for (i=0; i<npoints; i++) {
    Orig_pointset[i].point = (double*)malloc(dim*sizeof(double));
    Orig_pointset[i].index = i;
    for (j=0; j<dim; j++) {
      // newline counts as whitespace
      if (!fscanf(pointfile, "%lg ", &(Orig_pointset[i].point[j]))) {
	fprintf(stderr, "File does not contain enough data points!\n");
	exit(EXIT_FAILURE);
      }
    }
  }
  double super_best=1.0;
  for (nb_tries=0;nb_tries<10;nb_tries++){
  int rando1;
  int rando2;
  point swap;
  //SHUFFLING TIME: WORKS WITHOUT
  for (i=0;i<10*npoints;i++){
	  start=clock();
	  rando1= rand() % npoints;
	  rando2=rand() % npoints;
    swap = Orig_pointset[rando1];
    Orig_pointset[rando1] = Orig_pointset[rando2];
    Orig_pointset[rando2] = swap;
  }
  epsilon=1;
  for (i=0;i<npoints;i++){
	  for (j=i+1;j<npoints;j++){
		  if (fabs(Orig_pointset[i].point[0]-Orig_pointset[j].point[0])<epsilon)
			  epsilon=fabs(Orig_pointset[i].point[0]-Orig_pointset[j].point[0]);
	  }
  }
  best_bord=(double*)malloc(dim*sizeof(double));
  curr_bord=(double*)malloc(dim*sizeof(double));
  
  
	for (int i=0;i<dim;i++){ // Shouldn't be necessary as oydiscr necessarily modifies this at some point?
		best_bord[i]=0;
		curr_bord[i]=0;
	}
  optiset=(point*)malloc(kpoints*sizeof(point));  
  fill=false; // Should be useless
  for (i=0;i<kpoints;i++){
	  optiset[i].point=(double*)malloc(dim*sizeof(double));
	  memcpy(optiset[i].point,Orig_pointset[i].point, dim*sizeof(double));
    optiset[i].index = Orig_pointset[i].index;
  }
   
  upper=shift(npoints,kpoints,dim,Orig_pointset);
  end = clock();
  cput = ((double) (end - start)) / CLOCKS_PER_SEC;
  if (super_best>upper){
	  super_best=upper;
  FILE *fp; // Move our opti point set to a file
  fp=fopen(argv[5],"w");
  FILE *fpi;
  fpi=fopen(argv[6],"w");
  fprintf(fp, "n=%d,k=%d,dim=%d, discrepancy=%lf, runtime=%lf\n",npoints,kpoints,dim,upper,cput);
  for (i=0; i<kpoints;i++){
    fprintf(fpi,"%d\n",optiset[i].index);
	  for (j=0;j<dim;j++){
		  fprintf(fp,"%lf ",optiset[i].point[j]);
	  }
	  fprintf(fp,"\n");
  }
  fclose(fp);
  fclose(fpi);
  }
  for(i=0;i<kpoints;i++)
	  free(optiset[i].point);
  free(optiset);
  }
  fclose(pointfile);
  for (i=0;i<npoints;i++)
	  free(Orig_pointset[i].point);
  free(Orig_pointset);
  free(curr_bord);
  free(best_bord);
  fprintf(stderr,"We're done!");
  return 0;

}
