#ifndef UTILS_H
#define UTILS_H
/*****************************************************************************
 * 
 * File:    utils.h
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: February 2011
 *
 * Miscellaneous utilty functions
 *
 *
 ****************************************************************************/

#include <stdlib.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************
 *
 * Macros for debugging
 *
 ****************************************************************************/
  
#ifdef DEBUG_CONFIG
#define CONFIG_DEBUG_PRINT(x) printf("DEBUG CONFIG: "); printf x
#else
#define CONFIG_DEBUG_PRINT(x) /* nothing */
#endif
  
#ifdef DEBUG_DIGRAPH
#define DIGRAPH_DEBUG_PRINT(x) printf("DEBUG DIGRAPH: "); printf x
#else
#define DIGRAPH_DEBUG_PRINT(x) /* nothing */
#endif
  
#ifdef DEBUG_SAMPLER
#define SAMPLER_DEBUG_PRINT(x) printf("DEBUG SAMPLER: "); printf x
#else
#define SAMPLER_DEBUG_PRINT(x) /* nothing */
#endif

#ifdef DEBUG_ALGS
#define ALGS_DEBUG_PRINT(x) printf("DEBUG ALGS: "); printf x
#else
#define ALGS_DEBUG_PRINT(x) /* nothing */
#endif

#ifdef DEBUG_ALGEE
#define ALGEE_DEBUG_PRINT(x) printf("DEBUG ALGEE: "); printf x
#else
#define ALGEE_DEBUG_PRINT(x) /* nothing */
#endif

#ifdef DEBUG_MEMUSAGE
#define MEMUSAGE_DEBUG_PRINT(x) printf("DEBUG MEMUSAGE: "); printf x
#else
#define MEMUSAGE_DEBUG_PRINT(x) /* nothing */
#endif 
 

/*****************************************************************************
 *
 * Macros
 *
 ****************************************************************************/


#define FALSE 0
#define TRUE 1

#define MAX(a, b) ( (a) > (b) ? (a) : (b) )
  
/* Index into 2d n x n array stored row-major in contiguous memory  */
#define INDEX2D(i,j,n) ( ((i)*(n) + (j)) )

  
/*****************************************************************************
 *
 * typedefs
 *
 ****************************************************************************/

#ifdef UINT64_MAX
typedef unsigned long long uint64_t;
#endif
typedef unsigned int       uint32_t;
typedef unsigned int       uint_t;
typedef int                bool;
  

/*****************************************************************************
 *
 * function prototypes
 *
 ****************************************************************************/

  
/* pseudorandom numbers */

void init_prng(int tasknum); /* initialize the pseudorandom number generator */
double urand(void); /* uniform random double in [0,1] */
uint_t int_urand(uint_t n); /* uniform random int in 0...n-1 inclusive */
  
/* memory allocation */

void *safe_malloc(size_t size);
void *safe_calloc(size_t nelem, size_t elsize);
void *safe_realloc(void *ptr, size_t size);
char *safe_strdup(const char *s);

/* simple stats functions */

double mean_and_sd(double values[], uint_t nvalues, double *sd);

/* geographical functions */
  
long double deg2rad(long double deg);
long double rad2deg(long double rad);
double geo_distance(double lat1, double lon1, double lat2, double lon2);


/* miscellaneous */
  
int iDivUp(int a, int b);
int get_num_cores(void);
int timeval_subtract (struct timeval *result, struct timeval *x, 
                       struct timeval *y);
char *rstrip(char *s);
  
#ifdef __cplusplus
}
#endif

#endif /* UTILS_H */

