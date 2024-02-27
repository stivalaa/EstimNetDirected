/*****************************************************************************
 * 
 * File:    utils.c
 * Author:  Alex Stivala, Maksym Byshkin
 * Created: February 2011
 *
 * Miscellaneous utilty functions
 *
 *
 * Preprocessor defines used:
 *
 *    USE_POW_LOOKUP      - use lookup table for pow()
 *
 *
 ****************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#include "utils.h"

#ifdef USE_RANDOM123
#include <Random123/threefry.h>
#else
#error "Must use a good pseudorandom number generator"
#endif

/*****************************************************************************
 *
 * externally visible variables
 *
 ****************************************************************************/

#ifdef USE_POW_LOOKUP
#error "No longer in use as lambda is runtime parameter, and pow() is actually faster anyway"
/* Lookup table of integer power y of double x, faster than pow(x, y)
   intialized by init_powtable() */
double POWTABLE[POWTABLE_SIZE];
#endif

/*****************************************************************************
 *
 * pseudorandom numbers
 *
 ****************************************************************************/

#ifdef USE_RANDOM123
/* TODO don't use static variables so we can use threads */
static threefry2x64_ctr_t ctr = {{0, 0}};
static threefry2x64_key_t key = {{0xdeadbeef, 0xbadcafe}};
#endif

/*
 * Initialize the pseudorandom number genrator for given task number
 * (make sure seed is different for each task).
 */
void init_prng(int tasknum)
{
#ifdef USE_RANDOM123
  assert(sizeof(unsigned long long) == 8); /* since we use ULLONG_MAX */
  assert(sizeof(ctr.v[0]) == sizeof(unsigned long long));
  key.v[0] = time(NULL) + tasknum*123;
#else
#error "Must use a good pseudorandom number generator"
#endif
}

/*
 * Uniform random number in closed interval [0,1]
 */
double urand(void)
{
  /* TODO This is still not "really" uniform, although it is
     apparently what actual libraries use, due to non-uniform
     representation of reals by floating point, see
     http://mumble.net/~campbell/2014/04/28/uniform-random-float, also (older)
     better Doornik, J. A. (2007). Conversion of high-period random
     numbers to floating point. ACM Transactions on Modeling and
     Computer Simulation (TOMACS), 17(1), 3.
     http://www.doornik.com/research/randomdouble.pdf
     but it should be good enough */
#ifdef USE_RANDOM123
  double r;
  ctr.v[0]++;
  if (ctr.v[0] == 0) /* just in case we actually wrap 64 bit counter */
    ctr.v[1]++;
  threefry2x64_ctr_t randv = threefry2x64(ctr, key);
  r = (double)randv.v[0]/ULLONG_MAX;
  return r;
#else
#error "Must use a good pseudorandom number generator"
#endif
}

/*
 * Uniform random integer in 0..n-1 (inclusive)
 */
uint_t int_urand(uint_t n)
{
#ifdef USE_RANDOM123
  /* Make sure it is not biased
   * See https://petterhol.me/2020/09/01/common-mistakes-in-network-code/
   * https://www.pcg-random.org/posts/bounded-rands.html
   */
  ulonglong_t        threshold = -n % n;
  threefry2x64_ctr_t randv;
  ulonglong_t        r;
  do {
    ctr.v[0]++;
    if (ctr.v[0] == 0) /* just in case we actually wrap 64 bit counter */
      ctr.v[1]++;
    randv = threefry2x64(ctr, key);
    r =  randv.v[0];
  } while (r < threshold);
  return r % n;
#else
#error "Must use a good pseudorandom number generator"
#endif
}




/*****************************************************************************
 *
 * memory allocation
 *
 ****************************************************************************/


/*
 * malloc and exit on failure
 */
void *safe_malloc(size_t size)
{
  void *p = malloc(size);
  if (!p)  {
    fprintf(stderr, "malloc failed\n");
    exit(1);
  }
  return p;
}

void *safe_calloc(size_t nelem, size_t elsize)
{
  void *p = calloc(nelem, elsize);
  if (!p)  {
    fprintf(stderr, "calloc failed\n");
    exit(1);
  }
  return p;
}

void *safe_realloc(void *ptr, size_t size)
{
  void *p = realloc(ptr, size);
  if (size > 0 && !p)  {
    fprintf(stderr, "realloc() failed (%s)\n", strerror(errno));
    exit(1);
  }
  return p;
}

char *safe_strdup(const char *s)
{
  char *t;
  if (!(t = strdup(s))) {
    fprintf(stderr, "strdup failed\n");
    exit(1);
  }
  return t;
}


/*****************************************************************************
 *
 * simple stats functions
 *
 ****************************************************************************/


/*
 * Compute mean and standard deviation of an array of doubles
 *
 * Parameters:
 *    values  - array of values to compute mean and sd of
 *    nvalues - length of values array
 *    sd      - (Out) standard deviation
 *
 * Return value:
 *   mean of the values
 *
 */
double mean_and_sd(double values[], uint_t nvalues, double *sd)
{
  uint_t  i;
  double  mean = 0;
  double  meandiff;

  for (i = 0; i < nvalues; i++)
    mean += values[i];
  mean /= nvalues;
  *sd = 0;
  for (i = 0; i < nvalues; i++) {
    meandiff = values[i] - mean;
    *sd += meandiff * meandiff;
  }
  *sd = sqrt(*sd / nvalues);
  return mean;
}

/*****************************************************************************
 *
 * other utility functions
 *
 ****************************************************************************/


/*
 * iDivUp(a,b) = ceil(a / b) 
 */
int iDivUp(int a, int b)
{
  return ((a % b) != 0) ? (a / b + 1) : (a / b);
}


/*
 * Return the number of processors online
 */
int get_num_cores(void)
{
  long ncores;
  if ((ncores = sysconf(_SC_NPROCESSORS_ONLN)) < 0)  {
    fprintf(stderr,  "sysconf() failed (%d)\n", errno);
    exit(1);
  }
  return (int)ncores;
}


/* Subtract the `struct timeval' values X and Y,
   storing the result in RESULT.
   Return 1 if the difference is negative, otherwise 0.  
(from GNU libc manual) */
     
int
timeval_subtract (struct timeval *result, struct timeval *x, 
                  struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec) {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000) {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }
     
  /* Compute the time remaining to wait.
     tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;
     
  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

  
/* 
 * Strip whitespace from end of string in place
 */
char *rstrip(char *s)
{
  int i;
  for (i = strlen(s) - 1; i >= 0 && isspace(s[i]); i--)
    s[i] = '\0';
  return s;
}



/* compute three-dimensional Euclidean distance between two points with
   x,y,z coordinates */
double euclidean_distance(double x1, double y1, double z1,
                          double x2, double y2, double z2)
{
  return sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
}


/*
 * factorial(n)
 */
ulong_t factorial(ulong_t n)
{
  ulong_t fact = 1;
  while(n > 0) {
    fact *= n;
    fact--;
  }
  return fact;
}

/*
 * Initialize the integer power y of double x lookup table
 */
void init_powtable(double x)
{
#ifdef USE_POW_LOOKUP
  uint_t y;

  for (y = 0; y < POWTABLE_SIZE; y++) {
    POWTABLE[y] = pow(x, y);
  }
#endif
}



/*
 * Binomial coefficient n choose 2
 */
ulong_t n_choose_2(uint_t n)
{
  if (n < 2) {
    return 0;
  }
  return (ulong_t)n * (n - 1) / 2;
}

/*
 * Binomial coefficient n choose k
 */
ulonglong_t n_choose_k(uint_t n, uint_t k)
{
  uint_t i;
  double a = 1, b = 1;
  uint_t l = k;

  if (n < k) {
    return 0;
  }

  if (n - k < k) {
    l = n - k;
  }
  
  for (i = 1; i <= l; i++) {
    a *= (n + 1 - i);
    b *= i;
  }
  return a / b;
}


/*****************************************************************************
 *
 * geographical functions
 *
 ****************************************************************************/


const long double pi = 3.14159265358979323846;

/* convert degrees to radians */
long double deg2rad(long double deg)
{
  return deg * pi / 180;
}

/* convert radians to degrees */
long double rad2deg(long double rad)
{
  return rad * 180 / pi;
}

/* compute geographical (great-circle) distance in km betweeen two points
   specifed by latitude and longitude in degrees
   see e.g. https://en.wikipedia.org/wiki/Great-circle_distance*/
double geo_distance(double lat1, double lon1, double lat2, double lon2) {
  const double mean_earth_radius = 6371; /* km */
  long double theta, central_angle, dist;
  theta = lon1 - lon2;
  central_angle = acos( sin(deg2rad(lat1)) * sin(deg2rad(lat2))
                        + cos(deg2rad(lat1)) * cos(deg2rad(lat2))
                        * cos(deg2rad(theta)) );
  dist = mean_earth_radius * central_angle;
  return dist;
}

