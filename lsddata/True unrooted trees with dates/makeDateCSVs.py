import csv
import pandas as pd
import os

fns = [fn for fn in os.listdir( '.' ) if '.date' in fn ]

for datesfn in fns:
	x = pd.read_csv(datesfn,  header=0, delim_whitespace=True)
	dts = x.to_dict()[ x.to_dict().keys()[0] ]
	# pad names with zeros, len 6
	def pad_key( k ):
		l = len(k )
		return '0'*(6 - len(k)) + k
	#
	dts2 = dict( [ (pad_key(repr(k)), v) for k,v in dts.iteritems()] )
	
	f = open( datesfn + '.csv', 'w' )
	w= csv.writer( f )
	w.writerows( dts2.items() )
	f.close()
#
