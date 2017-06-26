'''
switch coalescent prior from constant to skyride 

elements to replace
1) beast -> 
	<!-- Generate a coalescent likelihood                                        -->
	<coalescentLikelihood id="coalescent">
		<model>
			<constantSize idref="constant"/>
		</model>
		<populationTree>
			<treeModel idref="treeModel"/>
		</populationTree>
	</coalescentLikelihood>

<!-- Generate a gmrfSkyrideLikelihood for GMRF Bayesian Skyride process      -->
<gmrfSkyrideLikelihood id="skyride" timeAwareSmoothing="true" randomizeTree="false">
	<populationSizes>

		<!-- skyride.logPopSize is in log units unlike other popSize                 -->
		<parameter id="skyride.logPopSize" dimension="SIZE" value="5.298317366548036"/>
	</populationSizes>
	<groupSizes>
		<parameter id="skyride.groupSize" dimension="SIZE"/>
	</groupSizes>
	<precisionParameter>
		<parameter id="skyride.precision" value="1.0" lower="0.0"/>
	</precisionParameter>
	<populationTree>
		<treeModel idref="treeModel"/>
	</populationTree>
</gmrfSkyrideLikelihood>

[ NOTE SIZE msut correspond to num taxa-1 !!] 




2) beast-> mcmc -> posterior -> prior -> 
				<oneOnXPrior>
					<parameter idref="constant.popSize"/>
				</oneOnXPrior>
<coalescentLikelihood idref="coalescent"/>

<gammaPrior shape="0.001" scale="1000.0" offset="0.0">
					<parameter idref="skyride.precision"/>
				</gammaPrior>
<gmrfSkyrideLikelihood idref="skyride"/>




3) beast -> operators -> 
		<scaleOperator scaleFactor="0.75" weight="3">
			<parameter idref="constant.popSize"/>
		</scaleOperator>
		
		<gmrfBlockUpdateOperator scaleFactor="2.0" weight="2">
			<gmrfSkyrideLikelihood idref="skyride"/>
		</gmrfBlockUpdateOperator>



4) beast -> mcmc -> log -> 
<parameter idref="constant.popSize"/>
<coalescentLikelihood idref="coalescent"/>

<parameter idref="skyride.precision"/>
			<parameter idref="skyride.logPopSize"/>
			<parameter idref="skyride.groupSize"/>
<gmrfSkyrideLikelihood idref="skyride"/>



5) beast -> mcmc
<mcmc id="mcmc" chainLength="10000000" autoOptimize="true">

<mcmc id="mcmc" chainLength="50000000" autoOptimize="true">
'''

import xml.etree.ElementTree as ET
import os

infns = os.listdir('beastxmls')

def make_skyride_lik( tree, size):
	SIZE = repr(size-1) 
	#~ gmrfSkyrideLikelihood  = ET.SubElement( tree, 'gmrfSkyrideLikelihood')
	gmrfSkyrideLikelihood  = ET.Element( 'gmrfSkyrideLikelihood')
	gmrfSkyrideLikelihood.set( 'id', 'skyride')
	gmrfSkyrideLikelihood.set('timeAwareSmoothing','true')
	gmrfSkyrideLikelihood.set('randomizeTree', 'false' )

	populationSizes = ET.SubElement( gmrfSkyrideLikelihood, 'populationSizes')
	parameter = ET.SubElement( populationSizes, 'parameter')
	parameter.set( 'id', 'skyride.logPopSize')
	parameter.set( 'dimension', SIZE)
	parameter.set( 'value', '4.0')

	groupSizes = ET.SubElement( gmrfSkyrideLikelihood, 'groupSizes')
	parameter1 = ET.SubElement( groupSizes, 'parameter')
	parameter1.set( 'id', 'skyride.groupSize')
	parameter1.set( 'dimension', SIZE )

	precisionParameter = ET.SubElement( gmrfSkyrideLikelihood, 'precisionParameter')
	parameter2 = ET.SubElement( precisionParameter, 'parameter')
	parameter2.set( 'id', 'skyride.precision')
	parameter2.set( 'value', '0.1')
	parameter2.set('lower', '0.0' )

	populationTree = ET.SubElement( gmrfSkyrideLikelihood, 'populationTree' )
	treeModel = ET.SubElement( populationTree, 'treeModel')
	treeModel.set( 'idref', 'treeModel')

	#~ return tree
	return gmrfSkyrideLikelihood
#


#~ infn = 'beastxmls/D750_3_25_62.xml'

for infn in infns:
	infn2 = 'beastxmls/' + infn 
	ofn = 'beastxmls.skyride/' + infn
	
	print infn2
	
	tree = ET.parse(infn2)
	beast = tree.getroot()

	size = len( beast.find('taxa').findall('taxon') )
	#~ skyridelik = make_skyride_lik( size ) 

	##1
	beast.remove( beast.find('coalescentLikelihood') )
	#~ beast = make_skyride_lik( beast, size)
	sk_lik_str  = ET.tostring(make_skyride_lik( beast, size)) # NOTE must insert this at proper place..

	##2
	prior = beast.find('mcmc').find('posterior').find('prior')
	prior.remove(  beast.find('mcmc').find('posterior').find('prior').find('oneOnXPrior') )
	prior.remove(  beast.find('mcmc').find('posterior').find('prior').find('coalescentLikelihood') )

	gammaPrior = ET.SubElement( prior, 'gammaPrior')
	gammaPrior.set( 'shape', '0.001')
	gammaPrior.set( 'scale', '1000.0')
	gammaPrior.set( 'offset', '0.0')
	parameter3 = ET.SubElement( gammaPrior, 'parameter')
	parameter3.set( 'idref', 'skyride.precision')
	gmrfSkyrideLikelihood = ET.SubElement( prior, 'gmrfSkyrideLikelihood')
	gmrfSkyrideLikelihood.set( 'idref', 'skyride')

	##3
	scaleops = beast.find('operators').findall('scaleOperator')
	scaleop0 = [ x for x in scaleops if x.find('parameter').get('idref')=='constant.popSize'][0]
	beast.find('operators').remove( scaleop0 )
	gmrfBlockUpdateOperator = ET.SubElement( beast.find('operators'), 'gmrfBlockUpdateOperator' )
	gmrfBlockUpdateOperator.set( 'scaleFactor', '2.0')
	gmrfBlockUpdateOperator.set( 'weight', '2')
	gmrfSkyrideLikelihood2 = ET.SubElement( gmrfBlockUpdateOperator, 'gmrfSkyrideLikelihood')
	gmrfSkyrideLikelihood2.set ( 'idref', 'skyride')


	##4
	#~ 4) beast -> mcmc -> log -> 
	#~ <parameter idref="constant.popSize"/>
	#~ <coalescentLikelihood idref="coalescent"/>
	filelog = beast.find('mcmc').findall('log')[1]
	logparms = filelog.findall('parameter')
	filelog.remove( [x for x in logparms if x.get('idref') == 'constant.popSize'][0] )
	filelog.remove( filelog.find('coalescentLikelihood') )

	parm5 = ET.SubElement( filelog, 'parameter')
	parm5.set( 'idref', 'skyride.precision')
	parm6 = ET.SubElement( filelog, 'parameter')
	parm6.set( 'idref', 'skyride.logPopSize')
	parm7 = ET.SubElement( filelog, 'parameter')
	parm7.set( 'idref', 'skyride.groupSize')
	gmrfSkyrideLikelihood4 = ET.SubElement( filelog, 'gmrfSkyrideLikelihood')
	gmrfSkyrideLikelihood4.set( 'idref', 'skyride')


	##5
	mcmc = beast.find('mcmc')
	mcmc.set('chainLength', '50000000' )

	##
	#~ tree.write('tmp3.xml' )

	##
	import re
	tree_str =  ET.tostring(tree.getroot()) 
	#~ </treeModel>
	#~ Definition:  re.sub(pattern, repl, string,
	tree_str2 = re.sub( '</treeModel>', '</treeModel>' + sk_lik_str, tree_str )


	f = open( ofn, 'w')
	f.write( tree_str2 )
	f.close()
