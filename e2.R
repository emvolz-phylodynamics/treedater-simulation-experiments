# experiment 2: make some tip dates uncertain
# see ../README
#R --args D750_11_10 F 11 T

PHI <- .20 # prop to be uncertain

require(treedater)

cargs <- commandArgs(trailingOnly=TRUE)
simcode <- cargs[1] #D750_11_10 D750_3_25 D995_11_10  D995_3_25
STRICT <- as.logical( cargs[2] )
WT <- as.numeric( cargs[3])
TEMPCONSTRAINT <- as.logical( cargs[4] ) 

OFN <- paste( sep='', 'e2/' , paste( sep='.', simcode, paste0('strict_',STRICT), WT
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

# mess up tip dates: 
{
lb <- range(sts)[1] 
ub <- range(sts)[2] 
nx <- round ( length( sts ) * PHI )
ix <- sample ( 1:length( sts ), size = nx, replace = F )
sts_lb <- sts 
sts_ub <- sts 
sts_lb[ix] <- lb 
sts_ub[ix] <- ub

est <- data.frame( lower = sts_lb, upper = sts_ub)
label_st_est <- names( sts )[ix] 
#~ est <- est[ label_st_est, ] 

stsd <- sd( sts ) 
.sts <- sts 
#~ .sts[ix] <- median( sts ) #NA
.sts[ix] <- pmin(ub, pmax(lb, rnorm( nx, sts[ix] , sd = stsd/2 ) ) )
}

print(paste(date(), 'start'))
runtime0 <- Sys.time()

td1 = td <- dater(tre, .sts[tre$tip.label], searchRoot=10, s=SEQLENGTH, omega0 = NA, maxit=300, temporalConstraints=TEMPCONSTRAINT, quiet=T,estimateSampleTimes = est )

runtime1 <- Sys.time()

# 2nd iteration
td2 <- dater(tre, td1$sts[tre$tip.label], searchRoot=10, s=SEQLENGTH, omega0 = NA, maxit=300, temporalConstraints=TEMPCONSTRAINT, quiet=T,estimateSampleTimes = est )

runtime2 <- Sys.time()

if (T)
{
print(paste(date(), 'dater', td$meanRate))
print(paste(date(), 'dater2', td2$meanRate))
pbtd = pbtd1 <- parboot.treedater( td1, nreps = 200 , overrideT=F)

runtime3 <- Sys.time()

pbtd2 <- parboot.treedater( td2, nreps = 200 , overrideT=F)

print(paste(date(), 'parboot'))
print( pbtd )

saveRDS( list(td = td1
	 , pbtd = pbtd1
	 , td2 = td2
	 , pbtd2 = pbtd2
	 , tre = tre 
	 , sts = sts 
	 , .sts = .sts
	 , label_st_est = label_st_est
	 , estimateSampleTimes = est 
	 , runtime0 = runtime0
	 , runtime1 = runtime1
	 , runtime2 = runtime2
	 , runtime3 = runtime3
 ), file = OFN 
)
}

