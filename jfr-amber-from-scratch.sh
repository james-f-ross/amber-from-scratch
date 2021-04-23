#!/bin/bash

# generate an amber simulation from scratch.
# only for cannonical amino acids.

# usage
# ./this.script input.pdb
# ./this.script 1PDB

# best if you supply a SSBOND pdb header, otherwise will attempt to autodetect 

########################################################################################
# INPUT OPTIONS
########################################################################################

# pdb can be a model present in the directory, if so, must include the .pdb extension
# pdb can just the pdb code, in which case it is fetched from the pdb database via pymol
pdb=$1

# choose simulation time (ns) # not scripted yet
time=1000

# location of amber/arc input templates
templatedir=~/work/8-scripts/1-MD-scripts/amber-from-scratch

# choose arc GPUs
arc=3p				# either 3k 3p or 4 		# 3k = arc3-K80, 3p = arc3-P100, 4 = arc4-V100

# define arc login and report details
USER=chmjro                 # string                # your username on the remote machine, if required


########################################################################################
# BODY TEXT - SETUP
########################################################################################

# check for user input
if [ -z $pdb ] ; then 
	echo "please provide a pdb ot pdb code"
	exit 0
fi

# set arc parameters
warc=$(echo $arc | cut -c 1)
if [ $warc = 3 ] ; then 
	COMPUTE=arc3.leeds.ac.uk
	kp=$(echo $arc | cut -c 2)
elif [ $warc = 4 ] ; then
	COMPUTE=arc4.leeds.ac.uk
	kp=''
else
	echo "arc compute undetermined"
	exit=0
fi 

# define arc report details and location 
EMAIL=$USER@leeds.ac.uk     # email address         # your email address
remotedir=/nobackup/$USER   # path                  # home directory on remote machine

# generate name variable
name=$(echo $pdb | sed 's/.pdb//g' )

# create working directory
oirdir=$(pwd) 
#mkdir amber-$name
#cd amber-$name
#cp ../$pdb . 2>/dev/null


# clean up the pdb in pymol
# removes everything that is not protein/nucleic acids
# removes alternate conformations
# removes anisotropic factor lines
if [ ! -f $pdb ] ; then 
	upname=$(echo $name | tr '[:lower:]' '[:upper:]')
	pdb=$upname.pdb
	name=$upname
	echo "fetch $upname, type=pdb1, async=0" > clean.pml
else
	echo "load $pdb" > clean.pml
fi 
echo "remove ! polymer
remove not alt ''+A
alter all, alt='' 
save $name.clean.pdb, $name, state=0" >> clean.pml
pymol -qc clean.pml
grep -v ANISOU $name.clean.pdb > $name.cleaner.pdb


# if its this big then we really ought to use high memory
highmem=no
if [ $(wc -l < $name.cleaner.pdb) -gt 200000 ] ; then 
	highmem=yes
fi

# if multistate structure, generate segments instead
count=1
rm $name.segment.pdb 2>/dev/null
for i in $(seq 1 $(grep MODEL $name.cleaner.pdb | wc -l ) ) ; do 
	MODEL=$(grep "MODEL.* $i$" $name.cleaner.pdb)
	if [ $i = 1 ] ; then  
	
		# check the number of segments
		nseg=$(sed -n "/^$MODEL/,/^ENDMDL/{p;/^ENDMDL/q}" $name.cleaner.pdb | grep ATOM | cut -c 73-76 | grep "\S" | sort | uniq | wc -l )
		if [ $nseg = 0 ] ; then 
			sed -i 's/./0/73;s/./0/74;s/./0/75;s/./0/76' $name.cleaner.pdb 
		fi
				
	fi
	for j in $(sed -n "/^$MODEL/,/^ENDMDL/{p;/^ENDMDL/q}" $name.cleaner.pdb | cut -c 73-76 | grep "\S" | sort | uniq ) ; do 
		segi=$(echo 000$count | rev | cut -c 1-4 | rev )
		sed -n "/^$MODEL/,/^ENDMDL/{p;/^ENDMDL/q}" $name.cleaner.pdb | grep "      $j " | sed "s/      $j /      $segi /g" >> $name.segment.pdb
		count=$(( $count + 1 ))
	done 
done


# add TER after each segment
count=1
cp $name.segment.pdb $name.TER.pdb
for i in $( cut -c 73-76 $name.segment.pdb | sort | uniq | grep "\S" ) ; do
for j in $( cut -c 22 $name.segment.pdb | sort | uniq | grep "\S" ) ; do
if [ $count = 1 ] ; then
count=$(( count + 1 ))
else 
count=$(( count + 1 ))
echo -en "\rAdding TER $count"
awk '/[A-Z][A-Z][A-Z] '$j' .* '$i' / && !x {print "TER"; x=1} 1' $name.TER.pdb > temp.pdb
mv temp.pdb $name.TER.pdb
fi
done
done
echo 'TER' >> $name.TER.pdb
uniq $name.TER.pdb > $name.TER2.pdb


# Looking for disulfides
# If downloaded, check pdb header.
native=$(echo $name | tr '[:upper:]' '[:lower:]')
if [ -f $native.pdb1 ] ; then 
	grep SSBOND $native.pdb1 > ssbond.tx
fi
echo ''

# If user suplied structure, check initial structure, otherwise measure sulphur distances 
if [ ! -f $native.pdb1 ] ; then 
	grep SSBOND $pdb > ssbond.tx
	if [ $(wc -l < ssbond.tx) = 0 ] ; then
		echo "auto detecting disulphide bonds"
		grep "SG  CY" $pdb | cut -c 18-27,31-38,39-46,47-54,72-76 | sed 's/CY/\nCY/g' | grep "\S" | sed 's/./&,/9;s/./&,/19;s/./&,/28;s/./&,/37' > sspos.tx 
		if [ $(wc -l < sspos.tx ) -gt 1 ] ; then
			awk -F',' '{print $1,$5}' sspos.tx | sed 's/ /_/g;s/__/_/g;s/__/_/g' > header.tx
			awk -F',' '{print $1,$5}' sspos.tx  > header.tx
			awk -F',' '{print $2,$3,$4}' sspos.tx | sed 's/ \+/,/g;s/\.//g;s/^,//g' > data.tx
			paste -d',' header.tx data.tx > hdata.tx
	
			echo "from sklearn.metrics.pairwise import euclidean_distances
import numpy as np 
X = np.loadtxt('hdata.tx',delimiter=',',skiprows=0,usecols=(1,2,3))
np.savetxt('distance-matrix2.csv', euclidean_distances(X,X),fmt='%i', delimiter=',') " > dmp.py
			python dmp.py
			wait
			paste -d',' header.tx distance-matrix2.csv > distance-matrix.csv

			sscount=0
			for i in $(seq 1 $(wc -l < header.tx) ) ; do
				count=$(( $i + 1 ))
				if [ $(awk -F',' -v var=$count '{print $1,$var}' distance-matrix.csv | grep " 2[0-9][0-9][0-9]$" | wc -l ) = 1 ] ; then 
					sscount=$(( $sscount + 1 ))
					ssnum=$(echo "   $sscount" | rev | cut -c 1-4 | rev)
					echo "SSBOND$ssnum $(sed -n ''$i'p' header.tx | cut -c 1-9)    $(awk -F',' -v var=$count '{print $1,$var}' distance-matrix.csv | grep " 2[0-9][0-9][0-9]$" | cut -c 1-9 )" | sed 's/./& /16;s/./& /30' >> ssbond.tx
					grep -v "$(sed -n ''$i'p' header.tx | cut -c 1-9)" distance-matrix.csv > temp.tx ; mv temp.tx distance-matrix.csv
					grep -v "$(awk -F',' -v var=$count '{print $1,$var}' distance-matrix.csv | grep " 2[0-9][0-9][0-9]$" | cut -c 1-9 )" distance-matrix.csv > temp.tx ; mv temp.tx distance-matrix.csv
				elif [ $(awk -F',' -v var=$count '{print $1,$var}' distance-matrix.csv | grep " 2[0-9][0-9][0-9]$" | wc -l ) -gt 1 ] ; then 
					echo "unable to determine disulphide bonds, please define mandually"
					exit 0 
				fi 
			done
		fi
	fi
fi

# report on disulphides 
if [ $(wc -l < ssbond.tx ) != 0 ] ; then 
	echo "found these SSBONDS"
	cat ssbond.tx
else	
	echo "did not find SSBONDS"
fi	
	
cat ssbond.tx $name.TER2.pdb > $name.ssbond.pdb

# pass to amber to check through structure.
# pdb4amber is broken for more than 9999 residues
if [ $( grep ' CA ' $name.ssbond.pdb | wc -l ) -lt 9999 ] ; then 
	pdb4amber -p $name.ssbond.pdb > $name.amber.pdb
else	
	cp $name.ssbond.pdb $name.amber.pdb
fi

########################################################################################
# BODY TEXT - EXPLICIT SOLVATON
########################################################################################

# set boundary for edge waters
gap=12.0

# set 

# estimate volume from min/max dimentions for salt concentration.
minx=$( grep ' CA ' $name.amber.pdb | cut -c 31-38 | sort -n | head -1 | sed 's/ //g' )
maxx=$( grep ' CA ' $name.amber.pdb | cut -c 31-38 | sort -n | tail -1 | sed 's/ //g' )
miny=$( grep ' CA ' $name.amber.pdb | cut -c 39-46 | sort -n | head -1 | sed 's/ //g' )
maxy=$( grep ' CA ' $name.amber.pdb | cut -c 39-46 | sort -n | tail -1 | sed 's/ //g' )
minz=$( grep ' CA ' $name.amber.pdb | cut -c 47-54 | sort -n | head -1 | sed 's/ //g' )
maxz=$( grep ' CA ' $name.amber.pdb | cut -c 47-54 | sort -n | tail -1 | sed 's/ //g' )
pvolume=$( echo "( ( $maxx + $gap ) - ( $minx - $gap ) ) * ( ( $maxy + $gap ) - ( $miny - $gap ) ) * ( ( $maxz + $gap ) - ( $minz - $gap ) ) * 0.666 " | bc | awk -F'.' '{print $1}')

# adjust volume by subtracting protein atoms - arbitratry factor of 5!
advolume=$(( $(tail $name.amber.pdb | awk '{print $2}' | grep "\S" | tail -1) * 10 ))
volume=$(( $pvolume - $advolume ))

# find 150 mM salt ions
ions=$(( $volume * 903 / 10000000 ))

# find charge to neutrilise
posi=$(grep ' CA ' $name.amber.pdb | grep "LYS\|ARG" | wc -l )
negi=$(grep ' CA ' $name.amber.pdb | grep "GLU\|ASP" | wc -l )
charge=$(( $posi - $negi ))
if [ $charge -gt 0 ] ; then
	salt="Cl-"
	nions=$(( $charge + $ions ))
	pions=$ions
elif [ $charge -lt 0 ] ; then
	salt="Na+"
	ncharge=$(( $negi - $posi ))
	pions=$(( $ncharge + $ions ))
	nions=$ions
else 
	pions=$ions
	nions=$ions
fi 
echo "volume=$volume therefore $ions ions : charge=$charge so Na=$pions and Cl=$nions"

# Generate solvated volume
rm leap.log 2>/dev/null
echo "source leaprc.protein.ff14SB
source leaprc.water.tip3p
mol = loadPdb $name.amber.pdb
solvateBox mol TIP3PBOX $gap
addIonsRand mol Na+ $pions Cl- $nions
savepdb mol $name.explicit.pdb
saveAmberParm mol $name.explicit.parm7 $name.explicit.rst7
quit" > tleap.solvate.in

# With the ff19SB force field there is an updated water model.
# can use : source leaprc.protein.ff19SB
# with    : source leaprc.water.opc
#		  : loadamberparams frcmod.ions1lm_126_hfe_opc
#		  : solvateOct ramp SPCBOX 14.0

tleap -f tleap.solvate.in

avolume=$(grep Volume leap.log | awk '{print $2}' | awk -F'.' '{print $1}')
mv leap.log leap.solvate.log
rm tleap.in 2>/dev/null

# compare calclutated volumes
echo "EstiVolume = $volume
CalcVolume = $avolume
Vol-Change = $(( $volume - $avolume ))"

# print volume charge log for future adjustment
echo "$name, dimensionV=$pvolume, pro-num=$advolume, EstiVolume = $volume, CalcVolume = $avolume, Vol-Change=$(( $volume - $avolume )), ions=$ions, charge=$charge so Na=$pions and Cl=$nions" >> $templatedir/salt-checker

echo "hmassrepartition
outparm $name.explicit.hmr.parm7
quit" > parm.in
parmed -i parm.in $name.explicit.parm7 

# generate unique names
rand=$RANDOM
date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')-$rand)

# make simulation folder and populate
mkdir $date-amber-$name

	# . . . with README's
cp $templatedir/README* $date-amber-$name/.

	# . . . with qsub submission scripts
if [ $warc = 3 ] ; then 
	sed "/^#$ -M.*/a #$ -N amin-$rand" $templatedir/run-amberA3-min.sh > $date-amber-$name/amin-$rand.sh
	if [ $highmem = yes ] ; then
		sed "/^#$ -M.*/a #$ -N amin-$rand" $templatedir/run-amberA3-mem-min.sh > $date-amber-$name/amin-$rand.sh
	fi
	sed "/^#$ -M.*/a #$ -N apro-$rand" $templatedir/run-amberA3-pro$kp.sh | sed "/^#$ -N apro-$rand.*/a #$ -hold_jid amin-$rand" > $date-amber-$name/apro-$rand.sh
	sed "/^#$ -M.*/a #$ -N apro-$rand" $templatedir/run-amberA3-rst$kp.sh | sed "/^#$ -N apro-$rand.*/a #$ -hold_jid amin-$rand" > $date-amber-$name/arst-$rand.sh.sh
elif [ $warc = 4 ] ; then
	sed "/^#$ -M.*/a #$ -N amin-$rand" $templatedir/run-amberA4-min.sh > $date-amber-$name/amin-$rand.sh
	if [ $highmem = yes ] ; then
		sed "/^#$ -M.*/a #$ -N amin-$rand" $templatedir/run-amberA4-mem-min.sh > $date-amber-$name/amin-$rand.sh
	fi
	sed "/^#$ -M.*/a #$ -N apro-$rand" $templatedir/run-amberA4-pro.sh | sed "/^#$ -N apro-$rand.*/a #$ -hold_jid amin-$rand" > $date-amber-$name/apro-$rand.sh
	sed "/^#$ -M.*/a #$ -N apro-$rand" $templatedir/run-amberA4-rst.sh | sed "/^#$ -N apro-$rand.*/a #$ -hold_jid amin-$rand" > $date-amber-$name/arst-$rand.sh.sh
else
	echo "arc compute undetermined"
	exit=0
fi 

	# . . . with md input files
sed "s/XXXX/$(grep ' CA ' $name.explicit.pdb | wc -l )/g" $templatedir/step4.0_minimization.mdin > $date-amber-$name/step4.0_minimization.mdin 
sed "s/XXXX/$(grep ' CA ' $name.explicit.pdb | wc -l )/g" $templatedir/step4.1_equilibration.mdin > $date-amber-$name/step4.1_equilibration.mdin 
cp $templatedir/step5_production-explicit.mdin $date-amber-$name/step5_production.mdin

	# . . . with structure files
mv $name.explicit.hmr.parm7 $date-amber-$name/step3_input.parm7
mv $name.explicit.rst7 $date-amber-$name/step3_input.rst7
mv $name.explicit.pdb $date-amber-$name/step3_input.pdb

	# pass to arc and submit jobs
ssh $USER@$COMPUTE mkdir $remotedir/$date-amber-$name
scp -r $date-amber-$name $USER@$COMPUTE:$remotedir/. 1>/dev/null
ssh $USER@$COMPUTE "cd $remotedir/$date-amber-$name ; qsub amin-$rand.sh ; qsub apro-$rand.sh" > $date.tx 


