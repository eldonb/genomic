
/* compile  - this one is for several homologous group of chromosomes - */
/*  need genomic coalescent simulator xibetagenomicpython.cpp */
/*  g++ -Wall -fPIC -shared -Ofast -o xibetapython.so xibetahomolmergerspython.cpp -lm -lgsl -lgslcblas */


using namespace std;
using namespace boost;
using namespace boost::random;
using namespace math;

/* obtain a seed out of thin air for the random number engine */
std::random_device randomseed;
  /* Standard mersenne twister  random number engine seeded with rng() */
  std::mt19937_64 rng(randomseed());
 /* construct a uniform on the unit interval; call with  runif( rng ) */
std::uniform_real_distribution<double> runif(0.0, 1.0);
std::uniform_int_distribution<unsigned> pchrom(1,4);
 /* generating a poisson random number generator */
std::poisson_distribution<unsigned> rpois(1.0) ;
std::exponential_distribution<double> rexp(1.0);
boost::random::beta_distribution<double> rbeta(1.0, 1.0); 
boost::random::uniform_smallint<unsigned> rsint(2,4); 
std::discrete_distribution<unsigned> rdist {};
std::binomial_distribution<unsigned> rbinom(1, .5); 

static double smallchoose( const double b, const unsigned k)
{

  double svar = 0 ;
  switch(k){
  case 0 : {svar = 1.; break; }
  case 1 : {svar = b; break; }
  case 2 : {svar = b*(b-1.)/2.; break;}
  case 3 : {svar =  b*(b-1.)*(b-2.)/6.; break;}
  case 4 : {svar = b*(b-1.)*(b-2.)*(b-3.)/24.; break;}
  default : break ;}
  return( svar );
  //return( (k < 1 ? 1. : (k < 2 ? b : ( k < 3 ? b*(b-1.)/2. : ( k < 4. ? b*(b-1.)*(b-2.)/6. : b*(b-1.)*(b-2.)*(b-3.)/24.) ) ) ) );
}


static double xreject( const double a, const double M )
{
  rbeta.reset();
  rbeta = boost::random::beta_distribution<double>(2.-a, a) ;
  double x;
  do{
    x = rbeta( rng);}
  while( (x > M) || (x <= 0.));
    
  assert( x > 0. );
  assert( x <=  M);
  return( x);
}

static double numerator(  std::forward_list<list<unsigned>>& l,   const double x)
{
  /* compute one term in the product over homologues */
  double y = 1. ;
  double n;
  for( auto &t: l){
    n = (double)t.size();
    /* y = y * (n > 1. ? ( pow(1. - x, n) +  (4.*n*(x/4.)*pow(1.-x, n-1.)) + (n > 1 ? 12.*smallchoose(n,2)*pow(x/4.,2)*pow(1.-x,n-2.) : 0.) + (n > 2 ?   24.*smallchoose(n,3)*pow(x/4.,3.)*pow(1.-x,n-3.) : 0.) + (n > 3 ? 24.*smallchoose(n,4)*pow(x/4.,4)*pow(1.-x,n-4.) : 0.)) : 1. );
     */
      y = y * (n > 1. ?  (pow(1. - x, n) +  (n*x*pow(1.-x, n-1.)) ) : 1.);

  }

  return( (1. - y)/(x*x) ); 
}


static unsigned sampletime( std::forward_list<list<unsigned>>& l , const double a, const double b,  const double M, std::vector<double>& w, double *Sjp,   double *xtimi )
{
  double rate,  exactf, z ;
  unsigned pwh; 
  rate = 0.;
  w.clear();
  for( auto& t: l){
    /* write n_l for number of lines at  lgroup l */
    z =  ( (unsigned)t.size() < 2 ? 0. : boost::math::binomial_coefficient<double>( (unsigned)t.size(), 2));
    w.push_back(z);
    /* rate is 4*choose( n_l, 2) */
    rate = rate + z; }
  /*   ratee = ratee + ((unsigned)t.size() < 3 ? 0. : boost::math::binomial_coefficient<double>( (unsigned)t.size(), 3));} */

  rdist.reset();
  rdist = std::discrete_distribution<unsigned> (w.begin(), w.end());
  /* draw time from  Exp(4*\sum_l choose(n_l, 2) ) */
  rexp.reset();
  rexp = std::exponential_distribution<double>(4.0*rate);
   
   do{
     xtimi[1] = xtimi[1] + rexp(rng); 
     xtimi[0] = (a < 2. ? xreject( a, M) : 0.);
     exactf = (xtimi[0] > 0.0 ? numerator( l, xtimi[0]) : rate);
     /* nm = (xtimi[0] < 0.0001 ? (rate - (2.*xtimi[0]*ratee)) : exactf)/rate; */
     /* nm =  exactf/rate; */
   }
   while( runif(rng) > (a < 2. ? exactf/rate : 0.25) );

   /* pick number of homolog seeing pairwise merger */
   pwh = rdist(rng);

    double Sj = (a < 2. ? 0. : (b > 0. ? log( exp(b*Sjp[pwh]) - (b*log(runif(rng)) / ( rate )))/b : Sjp[pwh] + xtimi[1] ) );
   xtimi[1] =  (a < 2. ? xtimi[1] : Sj - Sjp[pwh]) ;
   Sjp[pwh] = Sj ;

   
   return( pwh);
}

static unsigned samplebinom( const double n, const double x )
{
  rbinom.reset();
   rbinom = std::binomial_distribution<unsigned>(n, x);
   return( rbinom(rng) );
}


static unsigned  samplek( const double n, const double x)
{

  unsigned k ;
  if( n > 2.){
    rbinom.reset();
    rbinom = std::binomial_distribution<unsigned>(n-2., x);
    k = 2 + (x > 0.0 ? rbinom(rng) : 0);
    while( runif(rng) > 2./( (double)(k*(k-1))) ){
      k = 2 + (x > 0.0 ? rbinom(rng) : 0); }
    assert( k > 1);
    assert( k <= (unsigned)n);}
  else{
    k = 2 ; }

  return(k);
}



static bool inRange( const unsigned x, const unsigned l, const unsigned h)
{

  return(  (x - l) <=  (h-l));  
  /*  return( x >= Imin ? (x <= Imax ? 0 : ( x >= Jmin ? (x <= Jmax ? 1 : 3) : 3)) : 3);  */

}


/* update (X_I, X_J, S) statistic */
static void updatespectrumIJK( std::forward_list<vector<double>>& spectrum, std::forward_list<list<unsigned>>& tree, const double theta, const double timi )
{
  double m ;
  unsigned g; 
    if (theta > 0. ){
      rpois.reset();
      rpois = std::poisson_distribution<unsigned>(theta * timi);}

    std::forward_list<vector<double>>::iterator s; 
       s =  spectrum.begin();
     for( auto &t: tree){
       if( (unsigned)t.size() > 1){
	 /*  std::cout << (unsigned)t.size() << ' ' << activehoms << '\n'; */
	 for( auto &l: t){
	   /* sample new branch length or mutations */
	   m =  (theta > 0. ? rpois(rng) : timi) ;
	   /* add to estimate of total tree length or mutation */
	   g = (inRange(l, Kmin, Kmax) ? 2 : 3); 
	   (*s)[g] = (*s)[g] + m ;
	   /* add to spectrum */
	   g = (inRange( l, Imin, Imax) ? 0 : (inRange( l, Jmin, Jmax) ? 1 : 3) );
	   (*s)[ g ] =  (*s)[ g ] + m; }}
       s =  std::next(s, 1) ; }
}

static void prentatre(  std::forward_list<list<unsigned> >& tree  )
{
  /* print tree */
  std::cout << "------" << std::endl;
    for( auto &t: tree){
      for( auto &h: t){
	std::cout << h << ' ';}
      std::cout << std::endl ;}
  std::cout << "------" << std::endl;
}


static void prentaspectrum( const unsigned n,  std::forward_list<vector<double>> spectrum )
{

     for( auto& s: spectrum){
    /*     for( double& y: s){ std::cout << y  << ' ' ; } */
    for( unsigned j = 1 ; j < n ; j++){
      std::cout << s[j]/s[0] << '\n' ; }
    std::cout <<  '\n' ; }

}

/*
I think you might be looking at an old version of msprime. I remember
I had x / (4 - 4x) in a WIP commit at one stage, but as you say that's
wrong. The current implementation splits this into two parts: first we
sample a total participant number as 2 + Bin(n - 2, x) (plus a
rejection control step, line 5925), and then divide into four groups
with participation probability 1 / (4 - i) for the ith group, i = {0,
1, 2, 3} (line 5850). Beta-Xi-Sim is exactly the same in that regard.
*/


static bool null(  unsigned x )
{

  return( x > 0); 
}

static void untilmerger( const double a, const double M,   std::forward_list<list<unsigned> >& tree,  gsl_matrix_ushort * mcounts, std::vector<double>& w, double *xtimi )
{
  bool merger = 0;
  unsigned h = 0;
  unsigned pwh = 0;
  unsigned count = 0;
    /* update time until next merger and sample blocks until merger */
    while( merger == 0){
      /* update time, x from  Beta(2-a,a), and  pairwise lgroup */
      pwh = 1; /*   sampletime(tree, a, M, w,  xtimi ); */
      h = 0 ;
      for( auto &t: tree){
	if( (unsigned)t.size() > 1){
	/* at least two blocks; check if pairwise merger homolog */
	/* sample total number K of  blocks to try to  merge; if linkage group then  K = 2 + Binom(n-2, x), otherwise K=Binom(n-2, x)  */
	  count = (h == pwh ? samplek( (double)t.size(), xtimi[0]) : samplebinom( (double)t.size(), xtimi[0]) );
	  /* count1 = (h == pwh ? samplek( (double)t.size(), xtimi[0]) : 0); */
	  if( count > 1){
	    /* at least two blocks to assign uniformly with replacemtn  into four boxes */
	    /* mcounts  records the number in each box for each lgroup */
	    gsl_matrix_ushort_set( mcounts, h, 1,  samplebinom( (double)count,  0.25 ) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 1) > 1 ? 1 : 0));
	    /* count1 =  samplebinom( (double)count,  0.25 ); */
	    merger = ( merger ? merger : (gsl_matrix_ushort_get(mcounts, h, 1)  > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set(mcounts, h, 2, (count - gsl_matrix_ushort_get(mcounts, h, 1)  > 0 ? samplebinom( (double)(count - gsl_matrix_ushort_get(mcounts, h, 1)), 1./3.) : 0));
	    merger = ( merger ? merger : ( gsl_matrix_ushort_get(mcounts, h, 2) > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 2) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    gsl_matrix_ushort_set(mcounts, h, 3, count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get( mcounts, h, 2) > 0 ? samplebinom( (double)(count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get(mcounts, h,2)), 0.5 ) : 0);
	    merger = ( merger ? merger : ( gsl_matrix_ushort_get( mcounts, h,3) > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 3) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    gsl_matrix_ushort_set(mcounts, h, 4, count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get(mcounts, h,2) - gsl_matrix_ushort_get(mcounts,h,3));
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 4) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    merger = ( merger ? merger : ( gsl_matrix_ushort_get(mcounts,h,4) > 1 ? 1 : 0));}
	}
	h =h + 1;} }
}



static unsigned size(  std::list<unsigned>& t, const unsigned begin, const unsigned end )
{
  
  /*   std::cout << *(t.begin()) << ' ' <<  begin << ' ' << end  << std::endl; */
  
  assert( end - begin > 1 ); 
  std::list<unsigned>::iterator h_itr ;
 
  unsigned s= 0;
  for( h_itr = std::next( t.begin(), begin); h_itr !=  std::next( t.begin(), begin+end); ++h_itr){
    assert( *(h_itr) > 0);
    s = s + (*h_itr);
    /* label the block pointed to by hitr to be removed */
    *h_itr = 0; }
  assert( s > 1);
  return( s );
}

/* update tree and  return updated count of active lgroups */
static void updatetree(  std::forward_list<list<unsigned> >& tree,  gsl_matrix_ushort * mcounts, unsigned *activelg,  std::vector<unsigned>& v)
{
  unsigned  h = 0;
  unsigned size1 = 0;
  unsigned size2 = 0;
  unsigned size3 = 0;
  unsigned size4 = 0;
  
  list<unsigned>::iterator h_itr ;
  unsigned k = 0; 
  
     for( auto &t: tree){
       /* check if merging blocks  at lgroup h */
       if( gsl_matrix_ushort_get( mcounts, h, 0) > 0 ){
	    /* shuffle the  blocks at lgroup t  */
	    v.clear();
	    v.assign( t.begin(),  t.end() );
	    std::shuffle( v.begin(), v.end(), rng);
	    t.assign( v.begin(), v.end() );

	    /* check if merger in group 1 */
	    if ( gsl_matrix_ushort_get(mcounts,h,1) > 1){
	      size1 = 0;
	      /* size( t, 0, gsl_matrix_ushort_get(mcounts,h,1));  
	       */ 
	      for( h_itr = t.begin(); h_itr !=  std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1)); ++h_itr){
		assert( *h_itr > 0);
		size1 = size1 + (*h_itr);
		*h_itr = 0;}
		
	    }
	    else{ gsl_matrix_ushort_set( mcounts, h, 1, 0);}

	    
	    /* check if mergers in group 2 */ 
	    if ( gsl_matrix_ushort_get(mcounts,h,2) > 1){
	      size2 = 0 ;
	      /* size( t,  gsl_matrix_ushort_get(mcounts,h,1), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2));}
	       */
	      for( h_itr = std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1)); h_itr !=  std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2)); ++h_itr){
		assert( *h_itr > 0);
		size2 = size2 + (*h_itr);
		*h_itr = 0; }}
	    else{ gsl_matrix_ushort_set(mcounts, h, 2, 0);}

	    if( gsl_matrix_ushort_get( mcounts,h,3) > 1){
	      size3 = 0;
	      /*size( t,  gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2),   gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) + gsl_matrix_ushort_get(mcounts,h,3));} */
		for( h_itr = std::next(t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2)); h_itr != std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) + gsl_matrix_ushort_get(mcounts,h,3)); ++h_itr){
		  assert( *h_itr > 0);
		  size3 = size3 + (*h_itr);
		  *h_itr = 0;} }
	    else{ gsl_matrix_ushort_set(mcounts,h,3, 0) ;}

	    if( gsl_matrix_ushort_get(mcounts,h,4) > 1){
	      size4 = 0;
	      /*size( t,  gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3),  gsl_matrix_ushort_get(mcounts,h,1)  + gsl_matrix_ushort_get(mcounts,h,2) + gsl_matrix_ushort_get(mcounts,h,3) +  gsl_matrix_ushort_get(mcounts,h,4));} */
		  for( h_itr = std::next(t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3)); h_itr != std::next( t.begin(), gsl_matrix_ushort_get( mcounts,h,1) +  gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3) + gsl_matrix_ushort_get( mcounts,h,4)); ++h_itr){
		    k = k + 1;
		    assert( *h_itr > 0);
		    size4 = size4 + (*h_itr);
		    *h_itr = 0;}
		  assert( k ==  gsl_matrix_ushort_get(mcounts, h, 4) );}
	    else{ gsl_matrix_ushort_set(mcounts, h, 4, 0); }
	    
	    
	    /* add new  blocks */
	    if( gsl_matrix_ushort_get(mcounts,h,1) > 1){ assert( size1 > 1);  t.push_back(size1);}
	    if( gsl_matrix_ushort_get(mcounts,h,2) > 1){ assert( size2 > 1);  t.push_back(size2);}
	    if( gsl_matrix_ushort_get(mcounts,h,3) > 1){ assert( size3 > 1);  t.push_back(size3);}
	    if( gsl_matrix_ushort_get(mcounts,h,4) > 1){ assert( size4 > 1);  t.push_back(size4);}
	    /* remove blocks that merged */ 
	    t.remove(0);
	    /*  assert( std::any_of( t.begin(), t.end(), null) ); */
	    activelg[0] = activelg[0] - (t.size() > 1 ? 0 : 1);}
       h = h + 1; }
}


/* note the points are hardcoded - the statistic is (X_I, X_J, S) */
static void printsfs(   std::forward_list<vector<double>> sfs,  double *svar  )
{

  svar[0] = 0.0;
  svar[1] = 0.0;
  svar[2] = 0.0;

  
  for( const auto& s: sfs){ 
    assert( s[2] >= 0. ); 
    
    svar[0] = svar[0] +  (s[2] > 0.0 ? (s[0]/s[2]) : 0.0 ) ;
    svar[1] = svar[1] +  (s[2] > 0.0 ? (s[1]/s[2]) : 0.0 ) ;
    svar[2] = svar[2] + s[2];}

  svar[0] = svar[0]/( (double)lgroups );
  svar[1] = svar[1]/( (double)lgroups );
  svar[2] = svar[2]/( (double)lgroups );
}



/* check sampletime for b parameter beta for exponential growth */
extern "C" void xibetagenomic( const double a, const double theta, double *result)
{

  /* change into a > 2 if running kingman(beta) */
  const double b = 0.0 ;
  /*
   struct node{
    unsigned number;
    list<double> h;};
  */
  

  std::forward_list<list<unsigned> > tree;
  std::list<unsigned> homolog;
  std::list<unsigned>::iterator h_itr ;
  forward_list<list<unsigned>>::iterator tree_itr;
  forward_list<vector<double>> spectrum;
  forward_list<vector<double>>::iterator spectrum_itr;
  
  std::vector<unsigned> v ;
  std::vector<double> x, w;
  double * Sjp = (double *)calloc( lgroups,  sizeof(double)) ;
  

  /*
  const unsigned n = (unsigned)atoi(argv[1]);
  const unsigned H = (unsigned)atoi(argv[2]);
  */
  /* const double a = atof(argv[1]); */
  const double K = 0.0 ;
  /*   const double b = 0.0 ; */
  /*  const double theta = atof(argv[2]); */
  /*
    const unsigned trials = (unsigned)atoi(argv[3]);  */
  /* fix the seed */
  /*   rng.seed( (unsigned)atoi(argv[argc - 1]) );  */
  const double M = (a < 2. ? (K > 0. ? K/(K + 1. + (pow(2., 1. - a)/(a-1.))) : 1.) : 1. );

  /*
  std::cout << "sample size " << n << '\n';
  std::cout << "homologs " << H << '\n';
  std::cout << "alpha " << a << '\n';
  std::cout << "K " << K << '\n';
  std::cout << "M " << M << '\n';
 std::cout << "beta " << b << '\n';
 std::cout << "theta " << theta << '\n';
 std::cout << "trials " << trials << '\n';
 std::cout << "seed " << atoi(argv[8]) << '\n';
  */

  double * xtimi = (double *)calloc(2, sizeof(double));
  unsigned *activelg = (unsigned *)malloc( sizeof(unsigned) ); 
  gsl_matrix_ushort * mcounts = gsl_matrix_ushort_calloc( lgroups,  6); 

  unsigned merger = 0;
  unsigned h = 0;
  unsigned pwh = 0;
  unsigned count = 0;
  unsigned size1 = 0;
  unsigned size2 = 0;
  unsigned size3 = 0;
  unsigned size4 = 0;
  
  /* initialise spectrum */
  for( unsigned j = 0 ; j < lgroups ; j++){
    x.assign(samplesize, 0.);
    spectrum.push_front( x); }
  
  unsigned r = 0 ;
  while( r < trials ){
    r = r + 1 ;  
    /* initialise tree */
  tree.clear();
  for( unsigned j = 0 ; j < lgroups ; j++){
    homolog.assign(samplesize, 1);
    /* homolog.push_front( (double)n );  */
    tree.push_front( homolog); }
  
  activelg[0] = lgroups;
  xtimi[0] = 0.0;
  xtimi[1] = 0.0;
  while(activelg[0] > 0){
    /*  prentatre(tree); */
    /* sample time and x; initialise the weights for sampling the pairwise merger locus */
    /* untilmerger( a, M, tree, mcounts, w, xtimi ); */
    /* update time until next merger and sample blocks until merger */
    merger = 0 ;
    while( merger < 1){
      gsl_matrix_ushort_set_zero( mcounts); 
      /* update time, x from  Beta(2-a,a), and  pairwise lgroup */
      pwh =   sampletime(tree, a, b,  M, w, Sjp,   xtimi );
      h = 0 ;
      for( auto &t: tree){
	if( (unsigned)t.size() > 1){
	/* at least two blocks; check if pairwise merger homolog */
	/* sample total number K of  blocks to try to  merge; if linkage group then  K = 2 + Binom(n-2, x), otherwise K=Binom(n-2, x)  */
	  count = (h == pwh ? samplek( (double)t.size(), xtimi[0]) : (a < 2. ? samplebinom( (double)t.size(), xtimi[0]) : 0) );
	  /* assert( count <= (unsigned)t.size() ); */
	  /* count1 = (h == pwh ? samplek( (double)t.size(), xtimi[0]) : 0); */
	  if( count > 1){
	    /* at least two blocks to assign uniformly with replacemtn  into four boxes */
	    /* mcounts  records the number in each box for each lgroup */
	    gsl_matrix_ushort_set( mcounts, h, 1,  samplebinom( (double)count,  0.25 ) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 1) > 1 ? 1 : 0));
	    merger = ( merger > 0 ? merger : (gsl_matrix_ushort_get(mcounts, h, 1)  > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set(mcounts, h, 2, (count - gsl_matrix_ushort_get(mcounts, h, 1)  > 1 ? samplebinom( (double)(count - gsl_matrix_ushort_get(mcounts, h, 1)), 1./3.) : 0));
	    merger = ( merger > 0 ? merger : ( gsl_matrix_ushort_get(mcounts, h, 2) > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 2) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    gsl_matrix_ushort_set(mcounts, h, 3, (count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get( mcounts, h, 2) > 1 ? samplebinom( (double)(count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get(mcounts, h,2)), 0.5) : 0) );
	    /*  std::cout << 'c' << count << ' '  << gsl_matrix_ushort_get(mcounts, h,1) << ' ' << gsl_matrix_ushort_get(mcounts, h,2) << ' ' << gsl_matrix_ushort_get(mcounts, h,3) << std::endl;   */
	    /* assert( (unsigned)(gsl_matrix_ushort_get(mcounts, h,1) + gsl_matrix_ushort_get(mcounts, h,2) + gsl_matrix_ushort_get(mcounts, h,3)) <= (unsigned)t.size() );    */
	    merger = ( merger > 0 ? merger : ( gsl_matrix_ushort_get( mcounts, h,3) > 1 ? 1 : 0) );
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 3) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    gsl_matrix_ushort_set(mcounts, h, 4, count - gsl_matrix_ushort_get(mcounts, h,1) - gsl_matrix_ushort_get(mcounts, h,2) - gsl_matrix_ushort_get(mcounts,h,3));
	    gsl_matrix_ushort_set( mcounts, h, 0, (gsl_matrix_ushort_get( mcounts, h, 4) > 1 ? 1 : gsl_matrix_ushort_get( mcounts, h, 0) ));
	    merger = ( merger > 0 ? merger : ( gsl_matrix_ushort_get(mcounts,h,4) > 1 ? 1 : 0));} }
	h = h + 1; } }
    /*  sampled time and merger; update spectrum    */
    updatespectrumIJK( spectrum, tree, theta, xtimi[1] ); 
    xtimi[1] = 0.0 ;
    /* update tree  
    updatetree( tree, mcounts, activelg, v);
    */
    h = 0;
    for( auto &t: tree){
       /* check if merging blocks  at lgroup h */
       if( gsl_matrix_ushort_get( mcounts, h, 0) > 0 ){
	 /* assert( std::all_of( t.begin(), t.end(), null) ); */
	    /* shuffle the  blocks at lgroup t  */
	    v.clear();
	    v.assign( t.begin(),  t.end() );
	    std::shuffle( v.begin(), v.end(), rng);
	    t.assign( v.begin(), v.end() );
	    assert( std::all_of( t.begin(), t.end(), null) );

	    /* check if merger in group 1 */
	    if ( gsl_matrix_ushort_get(mcounts,h,1) > 1){
	      size1 = 0 ;
	      /* size( t, 0, gsl_matrix_ushort_get(mcounts,h,1));  
	       */
	      for( h_itr = t.begin(); h_itr !=  std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1)); ++h_itr){
		assert( *h_itr > 0);
		size1 = size1 + (*h_itr);
		*h_itr = 0;}
	      assert( size1 > 1);
	    }
	    else{ gsl_matrix_ushort_set( mcounts, h, 1, 0);}

	    
	    /* check if mergers in group 2 */ 
	    if ( gsl_matrix_ushort_get(mcounts,h,2) > 1){
	      size2 = 0;
	      /* size( t,  gsl_matrix_ushort_get(mcounts,h,1), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2));}
	       */
	      for( h_itr = std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1)); h_itr !=  std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2)); ++h_itr){
		assert( *h_itr > 0);
		size2 = size2 + (*h_itr);
		*h_itr = 0; }
	      assert( size2 > 1);
	    }	
	    else{ gsl_matrix_ushort_set(mcounts, h, 2, 0);}

	    if( gsl_matrix_ushort_get( mcounts,h,3) > 1){
	      size3 = 0; 
	      for( h_itr = std::next(t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2)); h_itr != std::next( t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) + gsl_matrix_ushort_get(mcounts,h,3)); ++h_itr){
		assert( *h_itr > 0);
		size3 = size3 + (*h_itr);
		*h_itr = 0;}
	      assert( size3 > 1);
	    }
	    else{ gsl_matrix_ushort_set(mcounts,h,3, 0) ;}

	    if( gsl_matrix_ushort_get(mcounts,h,4) > 1){
	      size4 = 0;  /*size( t,  gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3),  gsl_matrix_ushort_get(mcounts,h,1)  + gsl_matrix_ushort_get(mcounts,h,2) + gsl_matrix_ushort_get(mcounts,h,3) +  gsl_matrix_ushort_get(mcounts,h,4));}
			   */
		  for( h_itr = std::next(t.begin(), gsl_matrix_ushort_get(mcounts,h,1) + gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3)); h_itr != std::next( t.begin(), gsl_matrix_ushort_get( mcounts,h,1) +  gsl_matrix_ushort_get(mcounts,h,2) +  gsl_matrix_ushort_get(mcounts,h,3) + gsl_matrix_ushort_get( mcounts,h,4)); ++h_itr){
		    
		    assert( *h_itr > 0);
		    size4 = size4 + (*h_itr);
		    *h_itr = 0;}
		  assert( size4 > 1);
	    }
	    else{ gsl_matrix_ushort_set(mcounts, h, 4, 0); }	    
	    /* add new  blocks */
	    if( gsl_matrix_ushort_get(mcounts,h,1) > 1){ assert( size1 > 1);  t.push_back(size1);}
	    if( gsl_matrix_ushort_get(mcounts,h,2) > 1){ assert( size2 > 1);  t.push_back(size2);}
	    if( gsl_matrix_ushort_get(mcounts,h,3) > 1){ assert( size3 > 1);  t.push_back(size3);}
	    if( gsl_matrix_ushort_get(mcounts,h,4) > 1){ assert( size4 > 1);  t.push_back(size4);}
	    /* remove blocks that merged */ 
	    t.remove(0);
	    /* assert( std::all_of( t.begin(), t.end(), null) ); */
	    activelg[0] = activelg[0] - (t.size() > 1 ? 0 : 1);}
       h = h + 1;
    }
    /* finished updating tree */
    /* prentatre( tree) ;  */
  }
  /* all blocks at all lgroups have merged */
  printsfs( spectrum, result );

  }
  /* end of trials */

  
  /* free memory */
  free(Sjp);
  free( xtimi);
  gsl_matrix_ushort_free( mcounts);
  free( activelg);
  tree.clear(); 
  forward_list<list<unsigned> >().swap(tree);
  homolog.clear();
  list<unsigned>().swap( homolog);
  spectrum.clear();
  forward_list<vector<double>>().swap( spectrum);

  v.clear();
  std::vector<unsigned>().swap( v) ;
  x.clear();
  w.clear();
  std::vector<double>().swap(x);
  std::vector<double>().swap(w);
  
}
