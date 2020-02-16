
/* modules relevant for exponential and algebraic growth */

static void printmatrix( gsl_matrix *M, int radir, int dalkar )
{
  int i, j;
  for( i = 1; i <= radir; i++){
    for( j = 1; j <= dalkar; j++){
      printf( "%g ", gsl_matrix_get(M, i, j));}
    printf("\n");}
}


static void initmatrix( int leaves, int lgroups,  gsl_matrix_int *G )
{
  /* n is number of chroms per linkage group; h is number of linkage groups */
  int i,j;
  for( j = 1; j <= lgroups; j++){
    gsl_matrix_int_set(G,j, 0,  leaves);
    for( i = 1; i <= leaves ; i++){
      gsl_matrix_int_set(G, j,i, i);}}

 /* set the number of linkage groups  */
  gsl_matrix_int_set(G,0,0, lgroups );
}


static double expgrowthSj( int j,  double beta, double Sjp, gsl_rng *r)
{
  /* drawing time S_j for  exponential growth  */
  /* S_j = T_j + S_{j+1} where S_{n+1} = 0 and n sample size */
  
  assert( j > 1);
  double Sj = 0. ;
  /* make sure Sj becomes > Sjp so we get  positive interval time */
  do{ 
    Sj = (beta > 0. ? log( exp( beta * Sjp) -  (beta*log(gsl_rng_uniform(r))/gsl_sf_choose( j, 2)))/beta : Sjp + gsl_ran_exponential(r, 1./gsl_sf_choose( j, 2))) ;}
  while( Sj <= Sjp);

  assert( Sj > 0.);
  assert( Sj > Sjp);
  return( Sj);
}

static int expgrowthlinkagegroup( int n, int h, gsl_matrix_int *G, double beta, gsl_matrix *Stime,   gsl_rng *r )
{
  /* Stime is same  'S_{m+1}' time for all  trees; S_m is time when first see m-1 blocks */
  /* beta is exponential growth parameter */
  /* return the linkage group with smallest time */
  int j, ilargest, activeloci;
  double t, largest, Tm;
  /* make sure at least one linkage group still with at least two blocks */
  assert( gsl_matrix_int_get( G, 0,0) > 0 );
  /* initialize the new interval time stored in Stime[0,0] */
  gsl_matrix_set( Stime, 0,0, 0.);
  Tm = 0. ;
  do{ 
    do{
    ilargest = 1;
    largest = 0. ;
    activeloci = 0;
    j = 1;
  for( j = 1; j <= h; j++){
      assert( j <= h);
    if( gsl_matrix_int_get(G,j,0) > 1 ){
      activeloci = activeloci + 1;
      /* t = S_m for linkage group j */ 
      t =  expgrowthSj( gsl_matrix_int_get(G, j, 0), beta, gsl_matrix_get(Stime, j, gsl_matrix_int_get(G, j,0)+1), r);
      assert( t > 0. );
      Tm = t -  gsl_matrix_get(Stime, j, gsl_matrix_int_get(G, j,0)+1);
      /* compute T_m = S_{m+1} - S_m for linkage group j */
      /* check that t = S_j = T_j + S_{j+1} > S_{j+1} */
      assert( Tm > 0.);}
    else{
      t = 0. ;}
    ilargest = (j > 1 ? (t > 0. ? (largest < 1./Tm ? j : ilargest) : ilargest) : 1);
    largest = ( t > 0. ? (j == ilargest ? 1./Tm : largest) : largest);}
  /* sampled an interval time - add to total interval time */
  /* double check 
     assert( largest > 0. ); */
  assert( activeloci == gsl_matrix_int_get(G,0,0)); 
  assert( ilargest > 0);
  assert( ilargest <= h);
  assert( gsl_matrix_int_get( G, ilargest, 0) > 1);}
    while( largest <= 0. );
  gsl_matrix_set(Stime, 0,0, gsl_matrix_get(Stime, 0,0) + (1./largest));}
  /* add and draw new interval times until the two blocks merge */
  while( gsl_rng_uniform(r) > .25);
  /* return the index of  active linkage group with shortest time */
  /* ilargest is now the linkage group at which two blocks merge */
  /* record   the S_j time at locus ilargest; S_j is the time until first see j-1 blocks */
  /* assert( gsl_matrix_get( Stime, ilargest, gsl_matrix_int_get(G, ilargest, 0)) == 0. ); */
  gsl_matrix_set(Stime, ilargest,  gsl_matrix_int_get(G, ilargest, 0),  gsl_matrix_get( Stime, ilargest, gsl_matrix_int_get(G, ilargest, 0)+1) +  gsl_matrix_get(Stime,0,0));
  
  /* return the linkage group with shortest interval time */
  assert( ilargest <= h);
  return( ilargest);
}


static void picktwoblocks( int leaves,  gsl_matrix_int *G,  int lgroup, int *indexes,  int *blocks,   gsl_rng *r)
{
  /* lgroup is linkage group at which blocks merge */
  assert( gsl_matrix_int_get(G, lgroup, 0) > 1);
  int k = 1;
  int i = 0;

  blocks[0] = blocks[1] = 0;
  indexes[0] = indexes[1] = 0;
  /* collect all indexes of active blocks at linkage group lgroup */
  for( k = 1; k <= leaves; k++){
    if( gsl_matrix_int_get(G, lgroup, k) == k){
    indexes[i] = k;
    i = i + 1;}}

  assert( i > 1);
  assert( indexes[0] > 0);
  assert( indexes[1] > 0);

  /* sample two block indexes from  indexes of active blocks at lgroup */
  gsl_ran_choose(r, blocks, 2, indexes, gsl_matrix_int_get(G, lgroup, 0), sizeof(int));
  assert(blocks[0] > 0);
  assert(blocks[1] > 0);
}



static void growthmergeblocks( int ssize,  gsl_matrix_int *G, int pairwiselocus, int *pairblocks)
{
  /* n is sample size per linkage group */
  /* h is number of linkage groups */
  /* lgroup is linkage group for which to merge blocks */
  /* double-check if still segregating blocks at lgroup */
  int k;
  /* draw the linkage group seeing a pairwise merger */
  /* pairwiselocus = expgrowthlinkagegroup( samplesize, lgroups, ebeta,  &Stimi); */
  /* draw two blocks at lgroup picked for pairwise merger */
  /* picktwoblocks(G, pairwiselocus, indexes,   pairblocks, r); */
  /* assign the two blocks to one of four boxes; if in same box then merge blocks */
  /* the two blocks are in same box with prob 1/4 */
  /* if( gsl_rng_uniform(r) < .25 ) */
    /* the two blocks are in same box, so merge at locus pairwiselocus */
  /**/
    for( k = ( pairblocks[0] < pairblocks[1] ? pairblocks[0] : pairblocks[1]) ; k <= ssize; k++){
      gsl_matrix_int_set(G, pairwiselocus, k,  ((gsl_matrix_int_get( G, pairwiselocus, k) ==  pairblocks[0]) || (gsl_matrix_int_get( G, pairwiselocus, k) ==  pairblocks[1]) ?  (pairblocks[0] <  pairblocks[1] ? pairblocks[0] : pairblocks[1]) : gsl_matrix_int_get(G, pairwiselocus, k)));}
    /* update the count of active blocks at linkage group pairwiselocus */
    gsl_matrix_int_set(G, pairwiselocus, 0, gsl_matrix_int_get(G,pairwiselocus,0) - 1); 
    /* update number of  linkage groups with active  loci */
    gsl_matrix_int_set(G, 0, 0, gsl_matrix_int_get(G,0,0) - (gsl_matrix_int_get(G,pairwiselocus,0) == 1 ? 1 : 0));
  /* return if merger > 0 then at least one pairwise merger overall */
  /* run until at least one pairwise merger overall */
}

static void expgrowthgenealogy( int ssize, int lgroups,   double expbeta, gsl_matrix *Stimi,  gsl_matrix_int *G,  int *indexes,    int *bl,   gsl_matrix *BLS,   gsl_rng *r)
{
  /* run a genealogy */
  /* n is sample size, ie number of chromosomes per  linkage groups */
  /* h is number of linkage groups */
  /* G[0,0] is the number of groups with  segregating blocks, ie not reached mrca */

  initmatrix(ssize, lgroups, G);
  int pairwiselocus;
  gsl_matrix_set_zero(BLS);

  while( gsl_matrix_int_get(G,0,0) > 0){
    /* at least 1 linkage group with segregating blocks */ 
      /*static int expgrowthlinkagegroup( int n, int h, gsl_matrix_int *G, double beta, gsl_matrix *Stime,   gsl_rng *r ) */
      pairwiselocus =  expgrowthlinkagegroup( ssize, lgroups, G, expbeta, Stimi, r);
      picktwoblocks(ssize, G, pairwiselocus, indexes, bl, r);
    /* update BLS matrix of branch lengths after updating the forest */
      updatebls(ssize, lgroups, G, BLS, gsl_matrix_get(Stimi, 0,0));
    /* update the forest after the merger */
    /*  growthmergeblocks( int ssize,  gsl_matrix_int *G, int pairwiselocus, int *pairblocks, gsl_rng *r) */
    growthmergeblocks( ssize, G, pairwiselocus, bl);}
  printmatrix( Stimi, lgroups,  ssize+ 1);
}



static void runstatistics( int leaves, int lgroups, double expgrbeta,  double runs,  gsl_rng *r)
{
  /*run many genomealogies and record statistics */
  double b = 0.;
   gsl_matrix_int *G = gsl_matrix_int_calloc( lgroups + 1,  leaves + 1);
   gsl_matrix *Stimes = gsl_matrix_calloc( lgroups + 1,  leaves + 2);
   gsl_matrix *BLS = gsl_matrix_calloc( lgroups + 1,  leaves + 1);
   gsl_matrix *CBLS = gsl_matrix_calloc( lgroups + 1, leaves + 1);
    double * coords1 = (double *)calloc( lgroups + 1, sizeof( double));
  double * coords2 = (double *)calloc( lgroups + 2, sizeof( double));
   int *pairblocks = (int *)calloc( 2, sizeof(int));
  int *indexes = (int *)calloc( leaves, sizeof(int));
  double * coords = (double *)calloc( 2, sizeof( double));
  int j;
  while(b < runs){
    /**/
    /* static void expgrowthgenealogy( int ssize, int lgroups,   double expbeta, gsl_matrix *Stimi,  gsl_matrix_int *G,  int *indexes,    int *bl,   gsl_matrix *BLS,   gsl_rng *r) */
    gsl_matrix_set_zero( Stimes);
    expgrowthgenealogy( leaves, lgroups,  expgrbeta, Stimes, G,  indexes, pairblocks, BLS, r);
    /* sumbls(SAMPLESIZE, LINKAGEGROUPS, BLS, coords); */
    /*static void collectbls(int n, int h, gsl_matrix *BLS, gsl_matrix *collectedbls )*/
    /* collectbls(leaves, lgroups, BLS, CBLS); */
    /* static void sfsarray(int n, int H,  gsl_matrix * BLS, double *coord1, double *coord2) */
    sfsarray(leaves, lgroups, BLS, coords1, coords2);
     coords1[0] = 0.;
     coords2[0] = 0.;
     for( j = 1; j <=  lgroups ; j++){
       coords1[0] = coords1[0]  +  coords1[j] ;
       coords2[0] = coords2[0]  +  coords2[j] ;}
     printf("%g %g\n", coords1[0]/((double)lgroups),  coords2[0]/((double)lgroups));
    /* printf( "%g %g\n", coords[0], coords[1]); */
    
    b = b + 1.;}

  /* static void printm( int n, int h, double denom,  gsl_matrix *X)*/
  /* printm( leaves, lgroups, runs,  CBLS); */

  /* free mem */
   gsl_matrix_int_free( G);
   gsl_matrix_free(BLS);
   gsl_matrix_free(CBLS);
gsl_matrix_free(Stimes);
 free( coords1);
 free( coords2);
   free( pairblocks);
   free(indexes);
   free(coords);
}


static void printvbi( int leaves, int lgroups, double expgrbeta,  double runs,  gsl_rng *r)
{
  /*run many genomealogies and record statistics */
  double b = 0.;
   gsl_matrix_int *G = gsl_matrix_int_calloc( lgroups + 1,  leaves + 1);
   gsl_matrix *Stimes = gsl_matrix_calloc( lgroups + 1,  leaves + 2);
   gsl_matrix *BLS = gsl_matrix_calloc( lgroups + 1,  leaves + 1);
   gsl_matrix *CBLS = gsl_matrix_calloc( lgroups + 1, leaves + 1);
    double * coords1 = (double *)calloc( lgroups + 1, sizeof( double));
  double * coords2 = (double *)calloc( lgroups + 2, sizeof( double));
   int *pairblocks = (int *)calloc( 2, sizeof(int));
  int *indexes = (int *)calloc( leaves, sizeof(int));
  double * coords = (double *)calloc( 2, sizeof( double));
  int j, i;
  while(b < runs){
    /**/
    /* static void expgrowthgenealogy( int ssize, int lgroups,   double expbeta, gsl_matrix *Stimi,  gsl_matrix_int *G,  int *indexes,    int *bl,   gsl_matrix *BLS,   gsl_rng *r) */
    expgrowthgenealogy( leaves, lgroups,  expgrbeta, Stimes, G,  indexes, pairblocks, BLS, r);
    /* sumbls(SAMPLESIZE, LINKAGEGROUPS, BLS, coords); */
    /*static void collectbls(int n, int h, gsl_matrix *BLS, gsl_matrix *collectedbls )*/
    /* collectbls(leaves, lgroups, BLS, CBLS); */
    /* static void sfsarray(int n, int H,  gsl_matrix * BLS, double *coord1, double *coord2) */
    /* sfsarray(leaves, lgroups, BLS, coords1, coords2); */
    /*
     coords1[0] = 0.;
     coords2[0] = 0.;
    */
     for( j = 1; j <=  lgroups ; j++){
       for( i = 1; i < leaves; i++){
	 printf("%g ", gsl_matrix_get( BLS, j, i)); }
       printf( "\n"); }
     /*
       coords1[0] = coords1[0]  +  coords1[j] ;
       coords2[0] = coords2[0]  +  coords2[j] ;}
     printf("%g %g\n", coords1[0]/((double)lgroups),  coords2[0]/((double)lgroups));
      printf( "%g %g\n", coords[0], coords[1]);  */
    
    b = b + 1.;}
  /* static void printm( int n, int h, double denom,  gsl_matrix *X)*/
  /* printm( leaves, lgroups, runs,  CBLS); */

  /* free mem */
   gsl_matrix_int_free( G);
   gsl_matrix_free(BLS);
   gsl_matrix_free(CBLS);
gsl_matrix_free(Stimes);
 free( coords1);
 free( coords2);
   free( pairblocks);
   free(indexes);
   free(coords);
}
