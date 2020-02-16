poptreesizewf <-
    function( iN, iNe)
{
                                        # iN is pop size
                                        ## install.packages( pk="/home/BJARKI/verk/Rpakkar/Rcpp_1.0.2.tar.gz", repos=NULL, type="source") 
                                        ## install.packages( pk="/home/BJARKI/verk/Rpakkar/extraDistr", repos=NULL, type="source")
                                        ## library(extraDistr)

    numberanc <- iN
    
    size <- 0
    pop <- as.numeric(1:iNe)
                                        #current <- iN
    while(numberanc > 1){
        ##  yencoded
        size <- size + numberanc
        numberanc <- length(unique( sample( pop, size=numberanc,  replace=T ))) }
     
    return( size)
}
##################
"pow" <- function(x,y)
{
 return( (as.numeric(x))^(as.numeric(y)) )

}
###################
poptreesizeskewed <-
    function( iN, irangeX, iialpha,  skra)
{
                                        # iN is pop size
                                        ## install.packages( pk="/home/BJARKI/verk/Rpakkar/Rcpp_1.0.2.tar.gz", repos=NULL, type="source") 
                                        ## install.packages( pk="/home/BJARKI/verk/Rpakkar/extraDistr", repos=NULL, type="source")
                                        ## library(extraDistr)

    u <- iN
    lina <- 1
    size <- 0
                                        #current <- iN
    while(u > 1){
        size <- size + u
        ##  yencoded
        ## yencoded <-  sample( pop, size=u, replace=TRUE )
        ## sample the parent  levels
        ## first sample juveniles Xi for each of iN individuals
        
        iivpi <- pow( irangeX, -iialpha) - pow( irangeX +1, -iialpha)
        iivpi <- iivpi/sum(iivpi)
        vxi <- sample( irangeX, size=iN,  replace=T, prob=iivpi)
        ## shortcut: use the Xi in assigning lines to parents by multivariate hypergeom sampling
        lvui <- as.vector( rmvhyper( 1, vxi, u) )
        yencoded <- rep( 1:iN, lvui)
        ## print( yencoded)
        ## countanc <- table( parents)
        ##yencoded <- rep( 1:length(countanc), as.vector(countanc) )
                                        # parents  are not  in ascending order
                                        ##save(parents, file=paste(c("yencoded",lina,".rdata"),collapse=""), compression_level=1)
                                        ##  save(yencoded, file=paste(c(skra,lina,".rdata"),collapse=""), compression_level=1)
                                        #current <- unique(parents)
                                        #u <- length(current)
                                        #u <- length(unique(yencoded))
       ## print(yencoded)
        u <- length(unique(yencoded))
        lina <- lina + 1
    }
    return( size)
}


##############################
"teiknapoptreesize" <- function(itrials )
{

    sizes <- numeric(0)
    iiN <- 10^6
    ialpha <- 1.5

    safn <- numeric(0)
    for( i in 1:itrials){
        safn <- c( safn, poptreesizeskewed( iiN, as.numeric(1:(10^8)), ialpha, "dfdfd"))}

    sizes <- mean(safn)
    
    safn <- numeric(0)
   for( i in 1:itrials){
       safn <- c(safn,  poptreesizewf( iiN, pow( iiN, ialpha - 1) )) }

    sizes <- c( sizes, mean(safn))

    for( j in 1:6){
     safn <- numeric(0)
    for( i in 1:itrials){
        safn <- c(safn,  poptreesizewf( iiN, tpow(j) ) ) }
     sizes <- c( sizes, mean(safn) ) }

    print(sizes)
    
    print( sizes[2:8] / sizes[1])
    
   ## plot( sizes)
    
}
####################################
