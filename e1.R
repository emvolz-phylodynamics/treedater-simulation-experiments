# experiment 1: estimate rate, date, & parboot for CIs
# see ../README
# Rscript e1.R D750_11_10 T 15 T
# Rscript e1.R D750_11_10 F 15 T
# R --args D750_11_10 F 15 T

require(treedater)


cargs <- commandArgs(trailingOnly=TRUE)
simcode <- cargs[1] #D750_11_10 D750_3_25 D995_11_10  D995_3_25
STRICT <- as.logical( cargs[2] )
WT <- as.numeric( cargs[3])
TEMPCONSTRAINT <- as.logical( cargs[4] ) 

OFN <- paste( sep='', 'e1/' , paste( sep='.', simcode, paste0('strict_',STRICT), WT
  , paste0( 'tempConstraint_', TEMPCONSTRAINT )
  , 'rds') 
)

TRUTHDIR <- 'lsddata/True\ unrooted\ trees\ with\ dates/'

datadir <- ifelse( STRICT
 , 'lsddata/PhyML\ trees\ with\ BioNJ\ tree\ as\ strating\ tree/Strict\ molecular\ clock\ trees/'
 , 'lsddata/PhyML\ trees\ with\ BioNJ\ tree\ as\ strating\ tree/Relaxed\ molecular\ clock\ trees/'
)

SEQLENGTH <- 1e3
TRUEOMEGA <- 6e-3

#~ D750_11_10_rooted.tree    D750_3_25_rooted.tree    D995_11_10_rooted.tree    D995_3_25_rooted.tree
treefn <-  paste( datadir, simcode, '_unrooted.tree', sep='')
datesfn <- paste( TRUTHDIR, simcode, '.date' , sep='')


trees <- read.tree( treefn)

ststab <- read.table( datesfn , header=T)
sts <- setNames( ststab[,1], rownames( ststab))

tre <- unroot( trees[[WT]] )


print(paste(date(), 'start'))
runtime0 <- Sys.time()
td <- dater(tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=TEMPCONSTRAINT,  searchRoot = 10 , numStartConditions = 1)
#~ td <- dater(tre, sts[tre$tip.label], s=SEQLENGTH, omega0 = NA, maxit=300, quiet=T, temporalConstraints=TEMPCONSTRAINT,  searchRoot = 10 , numStartConditions = 0) #fast
runtime1 <- Sys.time()
print(paste(date(), 'dater', td$meanRate))
pbtd <- parboot.treedater( td, nreps = 200) # , override=T) 
#~ pbtd <- parboot.treedater( td, nreps = 100 , overrideT=T) #fast
runtime2 <- Sys.time()
print(paste(date(), 'parboot'))
print( pbtd )

#~ pbtd$trees <- NULL # comment this out to reduce storeage reqd

if (T)
{
saveRDS( list(td = td
	 , pbtd = pbtd
	 , tre = tre 
	 , sts = sts 
	 , runtime0 = runtime0
	 , runtime1 = runtime1
	 , runtime2 = runtime2
 ), file = OFN  # 
)
}



