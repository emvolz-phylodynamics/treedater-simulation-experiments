# check relaxed clock test
# Rscript e1.R D750_11_10 T 15 T
# R --args D750_11_10 F 15 T

require(treedater)


cargs <- commandArgs(trailingOnly=TRUE)

STRICT <- as.logical( cargs[1] )


SCENS <- c('D750_11_10', 'D750_3_25', 'D995_11_10' , 'D995_3_25' )


TRUTHDIR <- 'lsddata/True\ unrooted\ trees\ with\ dates/'

datadir <- ifelse( STRICT
 , 'lsddata/PhyML\ trees\ with\ BioNJ\ tree\ as\ strating\ tree/Strict\ molecular\ clock\ trees/'
 , 'lsddata/PhyML\ trees\ with\ BioNJ\ tree\ as\ strating\ tree/Relaxed\ molecular\ clock\ trees/'
)

SEQLENGTH <- 1e3
TRUEOMEGA <- 6e-3

#~ D750_11_10_rooted.tree    D750_3_25_rooted.tree    D995_11_10_rooted.tree    D995_3_25_rooted.tree
if (F)
{
lapply( SCENS, function( simcode ){
	treefn <-  paste( datadir, simcode, '_unrooted.tree', sep='')
	datesfn <- paste( TRUTHDIR, simcode, '.date' , sep='')
	
	trees <- read.tree( treefn)
	
	ststab <- read.table( datesfn , header=T)
	sts <- setNames( ststab[,1], rownames( ststab))
	
	sapply( trees, function(tre){
		rct <- relaxed.clock.test( tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=FALSE,  searchRoot = 5 , numStartConditions = 0)
		rct$clock
	})
} ) -> results

saveRDS( results, paste0('e3-strict', STRICT, '.rds') )
}

#~ STRICT == TRUE
#~ > summary( as.numeric(  unlist( results ) == 'strict' )  ) 
#~    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#~  0.0000  0.0000  1.0000  0.6525  1.0000  1.0000 

#~ STRICT==FALSE
#~ > summary( as.numeric(  unlist( results ) == 'strict' )  ) 
#~    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#~       0       0       0       0       0       0 


# inspect scen 2 
if (F)
{
simcode <- SCENS[2] 
	treefn <-  paste( datadir, simcode, '_unrooted.tree', sep='')
	datesfn <- paste( TRUTHDIR, simcode, '.date' , sep='')
	
	trees <- read.tree( treefn)
	
	ststab <- read.table( datesfn , header=T)
	sts <- setNames( ststab[,1], rownames( ststab))
	
	tre <- trees[[1]]
	rct <- relaxed.clock.test( tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=FALSE,  searchRoot = 5 , numStartConditions = 0)
	
	rct <- relaxed.clock.test( tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=TRUE,  searchRoot = 5 , numStartConditions = 0)
	
	rct$clock
}

# char all tree
if (F)
{
bdts <- lapply( SCENS, function(simcode ){
	treefn <-  paste( datadir, simcode, '_unrooted.tree', sep='')
	datesfn <- paste( TRUTHDIR, simcode, '.date' , sep='')
	
	trees <- read.tree( treefn)
	
	ststab <- read.table( datesfn , header=T)
	sts <- setNames( ststab[,1], rownames( ststab))
	
	tre <- trees[[1]]
	dater( tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=FALSE,  searchRoot = 5 , numStartConditions = 0)
})

#~ for ( bdt in bdts ){
for ( k in 1:4){
	X11(); plot ( ladderize(bdts[[k]]), main = SCENS[k] )
}
}
