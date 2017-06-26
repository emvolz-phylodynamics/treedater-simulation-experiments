require(ape)

infns <- list.files( pattern='.tree', rec = TRUE , full.names=T)
for (fn in infns )
{
	tres <- read.tree( fn )
	for (i in 1:length(tres )){
		write.tree( tres[[i]], file = paste0( fn, '_', i-1, '.nwk'))
	}
}
