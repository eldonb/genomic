#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mcheck.h>
/* 
#include <gmp.h> 
#include <mpfr.h>
#include <mpc.h>
*/
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_elementary.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_cdf.h>
#include </home/BJARKI/verk/truncation/wfiles/statistics.h>
#include </home/BJARKI/verk/truncation/wfiles/growthmodules.h>


static const double EXPGROWTH = 0.1;
static const double NUMBRUNS = 10.0;

enum{
  SAMPLESIZE = 10};

enum{
  LINKAGEGROUPS = 4};

gsl_rng *rngtype;

static void setup_rng(unsigned long int s)
{
const gsl_rng_type *T;
gsl_rng_env_setup();
T= gsl_rng_default;
rngtype= gsl_rng_alloc(T);
gsl_rng_set(rngtype,s);
}


int main(int argc, char *argv[])
{
  setup_rng(atol(argv[argc-1]) );
  /* static void runstatistics( int leaves, int lgroups, double expgrbeta,  double runs,  gsl_rng *r)  */
     runstatistics( SAMPLESIZE, LINKAGEGROUPS, EXPGROWTH, NUMBRUNS, rngtype ); 
  /* static void printvbi( int leaves, int lgroups, double expgrbeta,  double runs,  gsl_rng *r) */
  /* printvbi( SAMPLESIZE, LINKAGEGROUPS, EXPGROWTH,  NUMBRUNS, rngtype); */

  gsl_rng_free(rngtype);

  return GSL_SUCCESS;
}
