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
# pdb can just the pdb code, in which case it is fetched from the pdb database via /mnt/c/Users/james/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe
pdb=$1

# choose simulation time (ns) # not scripted yet
time=1000

# Salt concentration (M)
saltconc=0.150

# choose arc GPUs
arc=4                        # either 3k 3p or 4             # 3k = arc3-K80, 3p = arc3-P100, 4 = arc4-V100

# define arc login and report details
USER=chmjro                 # string                # your username on the remote machine, if required

# are you providing a restaint file? if so, name it
rst=rst					# filename or False

########################################################################################
# Generate folders and names
########################################################################################

# create working directory
oridir2=$(pwd)
oridir=$(echo $(pwd) | sed 's/ /\\ /g')
# generate name variable
name=$(echo $pdb | sed 's/.pdb//g' )
# generate unique names
rand=$RANDOM
date=$(echo $(date -R | cut -c 6-7,9-11,15-25 | sed 's/ /-/g;s/://g')-$rand)
# make simulation folder
outdir=$oridir/$date-amber-$name
mkdir $outdir

# define arc report details and location
EMAIL=$USER@leeds.ac.uk     # email address         # your email address
remotedir=/nobackup/$USER   # path                  # home directory on remote machine
PROXYJUMP=chmjro@remote-access.leeds.ac.uk  

# prepare proxyjump
if [ $PROXYJUMP = FALSE ] ; then
	PROXY=''
else
	PROXY=-oProxyJump=$PROXYJUMP
fi

########################################################################################
# DECLARE NON-DEFAULT INPUT OPTIONS
########################################################################################

# declare number of parts # dont change these please.
nmin=1
neqi=1
npro=1

# EG. For second equilibration step with smaller restaints, state additional parts above the first and use a '_neqi' suffix to specify changed parameters
#       neqi=2
#       restraint_wt_neqi2=1.0

# EG. for queuing 10 additional prodution runs (as default each run is 100 ns, so below generates 1 us, (provided 100 ns is attained in 48 hours)
# note the use of the '_npro' suffix in the variable, this is necessary.
npro=10
nstlim_npro2=50000000
nstlim_npro3=75000000
nstlim_npro4=100000000
nstlim_npro5=125000000
nstlim_npro6=150000000
nstlim_npro7=175000000
nstlim_npro8=200000000
nstlim_npro9=225000000
nstlim_npro10=250000000




########################################################################################
# INPUT SIMULATION OPTIONS (DEFAULTS) only change if you want to change the defaults
########################################################################################

# Restraint mask examples
# ='@CA,C,O,N&!:WAT' 	# backbone atoms and not water
# =':1-603'				# atom?/residue? range

# can use rst in cpptraj to generate harmonic distance restraints for nmropt options
# https://amberhub.chpc.utah.edu/rst/
# print out a default minimisation input
# Minimization input file in explicit solvent

# Defualt Minimization options (try not to change!)
imin=1         # Turn on minimization
maxcyc=5000    # Maximum number of minimization cycles
ncyc=2500      # 100 steepest-descent steps, better for strained systems
# Potential energy function options
cut=12.0       # nonbonded cutoff, in Angstroms
fswitch=10.0   # Force-based switching
# Control how often information is printed to the output file
ntpr=100       # Print energies every 100 steps
ntxo=2         # Write NetCDF format
# Restraint options
ntr=1          # Positional 07Jan22-140951-29497-amber-4crz
#restraints for proteins, sugars, and ligands
restraint_wt=1.0
restraintmask='@CA,C,O,N&!:WAT'
# Set water atom/residue names for SETTLE recognition
watnam='WAT'    # Water residues are named WAT
owtnm='O'       # Water oxygens are named O

for i in $(seq 1 $nmin) ; do
# redefine interative variables
eval  "$( ( set -o posix ; set ) | grep nmin$i= | sed 's/_nmin'$i'//g') "

# print mdin file
echo "Minimisation run
 &cntrl
    imin=$imin, maxcyc=$maxcyc, ncyc=$ncyc,
    cut=$cut, fswitch=$fswitch,
    ntpr=$ntpr, ntxo=$ntxo,
    ntr=$ntr, restraint_wt=$restraint_wt, restraintmask='$restraintmask',
    watnam='$watnam', owtnm='$owtnm',
 /
 &ewald
    vdwmeth = 0,
 /
END" > min$i.mdin
done

# print out a default equilibration input
# A NVT simulation for common production-level simulations

# Defualt equilibration options (try not to change!)
imin=0        # No minimization
irest=0       # This is NOT a restart of an old MD simulation
ntx=1         # So our inpcrd file has no velocities
# Temperature control
ntt=3         # Langevin dynamics
gamma_ln=1.0  # Friction coefficient (ps^-1)
tempi=303.15   # Initial temp -- give it some small random velocities
temp0=303.15   # Target temperature
# Potential energy control
cut=12.0      # nonbonded cutoff, in Angstroms
fswitch=10.0  # Force-based switching
# MD settings
nstlim=125000 # 125K steps, 125 ps total
dt=0.001      # time step (ps)
# SHAKE
ntc=2         # Constrain bonds containing hydrogen
ntf=2         # Do not calculate forces of bonds containing hydrogen
# Control how often information is printed
ntpr=1000     # Print energies every 1000 steps
ntwx=5000     # Print coordinates every 5000 steps to the trajectory
ntwr=10000    # Print a restart file every 10K steps (can be less frequent)
ntxo=2        # Write NetCrestraint_wt=1.0DF format
ioutfm=1      # Write NetCDF format (always do this#)
# Wrap coordinates when printing them to the same unit cell
iwrap=0
# Restraint options
ntr=1         # Positional restraints for proteins, sugars, and ligands
restraint_wt=1.0
restraintmask='@CA,C,O,N&!:WAT'
# Set water atom/residue names for SETTLE recognition
watnam='WAT'  # Water residues are named WAT
owtnm='O'     # Water oxygens are named O

for i in $(seq 1 $neqi) ; do
# redefine interative variables
eval  "$( ( set -o posix ; set ) | grep neqi$i= | sed 's/_neqi'$i'//g') "

# print mdin file
echo "Equilibration run
 &cntrl
    imin=$imin, irest=$irest, ntx=$ntx,
    ntt=$ntt, gamma_ln=$gamma_ln, tempi=$tempi, temp0=$temp0,
    cut=$cut, fswitch=$fswitch,
    nstlim=$nstlim, dt=$dt, ntc=$ntc, ntf=$ntf,
    ntpr=$ntpr, ntwx=$ntwx, ntwr=$ntwr, ntxo=$ntxo, ioutfm=$ioutfm, iwrap=$iwrap,
    ntr=$ntr, restraint_wt=$restraint_wt, restraintmask='$restraintmask',
    watnam='$watnam', owtnm='$owtnm',
 /
 &ewald
    vdwmeth = 0,
 /
END" > eqi$i.mdin

done


# print out a default production input

# Defualt production options (try not to change!)
imin=0        # No minimization
irest=1       # This IS a restart of an old MD simulation
ntx=5        # So our inpcrd file has velocities
if [ rst != False ] ; then 
nmropt=1	# restaints file included?? 1 yes, 0 no
else
nmropt=0
fi 
# Temperature control
ntt=3        # Langevin dynamics
gamma_ln=1.0  # Friction coefficient (ps^-1)
temp0=303.15   # Target temperature
# Potential energy control
cut=12.0      # nonbonded cutoff, in Angstroms
fswitch=10.0  # Force-based switching
# MD settings
nstlim=25000000 # 100 ns total
dt=0.004     # time step (ps)
ntc=2         # Constrain bonds containing hydrogen
ntf=2         # Do not calculate forces of bonds containing hydrogen
# Control how often information is printed
ntpr=25000    # Print energies every 100 ps
ntwx=125000   # Print coordinates every 500 ps to the trajectory
ntwr=125000   # Print a restart file every 500ps
ntxo=2        # Write NetCDF format
ioutfm=1      # Write NetCDF format (always do this#)
iwrap=1       # Wrap coordinates when printing them to the same unit cell
# Restraint options
ntr=0         # Positional restraints for proteins, sugars, and ligands
restraint_wt=0.0
restraintmask='@CA,C,O,N&!:WAT'
# Constant pressure control.
barostat=2    # MC barostat... change to 1 for Berendsen
ntp=1         # 1=isotropic, 2=anisotropic, 3=semi-isotropic w/ surften
pres0=1.0     # Target external pressure, in bar
# Set water atom/residue names for SETTLE recognition
watnam='WAT'  # Water residues are named WAT
owtnm='O'     # Water oxygens are named O

for i in $(seq 1 $npro) ; do
# redefine interative variables
eval  "$( ( set -o posix ; set ) | grep npro$i= | sed 's/_npro'$i'//g') "

# print mdin file
echo "Production run
 &cntrl
    imin=$imin, irest=$irest, ntx=$ntx, nmropt=$nmropt
    ntt=$ntt, gamma_ln=$gamma_ln, temp0=$temp0,
    cut=$cut, fswitch=$fswitch,
    nstlim=$nstlim, dt=$dt, ntc=$ntc, ntf=$ntf,
    ntpr=$ntpr, ntwx=$ntwx, ntwr=$ntwr, ntxo=$ntxo, ioutfm=$ioutfm, iwrap=$iwrap,
    ntr=$ntr, restraint_wt=$restraint_wt, restraintmask='$restraintmask',
    barostat=$barostat, ntp=$ntp, pres0=$pres0,
    watnam='$watnam', owtnm='$owtnm',
 /
 &ewald
    vdwmeth = 0,
 /
" > pro$i.mdin
if [ $nmropt = 1 ] ; then 
echo "DISANG=$rst" >> pro$i.mdin
fi
echo "END" >> pro$i.mdin

done

########################################################################################
# INPUT SIMULATION SUBMISSIONS (DEFAULTS)
########################################################################################

# submission scripts, minimisations
for i in $(seq 1 $nmin) ; do
echo "#!/bin/csh
mpirun pmemd.MPI -O -i min$i.mdin -p initial.parm7 -c initial.rst7 -o min$i.mdout -r min$i.rst7 -inf min$i.mdinfo -ref initial.rst7
end" > minq$i
pre=min$i
done

# submission scripts, equilibrations
for i in $(seq 1 $neqi) ; do
echo "#!/bin/csh
pmemd.cuda_SPFP -O -i eqi$i.mdin -p initial.parm7 -c $pre.rst7 -o eqi$i.mdout -r eqi$i.rst7 -inf eqi$i.mdinfo -ref initial.rst7 -x eqi$i.nc
end" > eqiq$i
pre=eqi$i
done

# submission scripts, productions
for i in $(seq 1 $npro) ; do
echo "#!/bin/csh
pmemd.cuda_SPFP -O -i pro$i.mdin -p initial.parm7 -c $pre.rst7 -o pro$i.mdout -r pro$i.rst7 -inf pro$i.mdinfo -ref initial.rst7 -x pro$i.nc
end" > proq$i
pre=pro$i
done



########################################################################################
# BODY TEXT - SETUP
########################################################################################

# check for user input
if [ -z $pdb ] ; then
        echo "please provide a pdb or pdb code"
        exit 0
fi

# set arc parameters
warc=$(echo $arc | cut -c 1)
if [ $warc = 3 ] ; then
        COMPUTE=arc3.leeds.ac.uk
        kp=$(echo $arc | cut -c 2)
elif [ $warc = 4 ] ; then
        COMPUTE=arc4.leeds.ac.uk
        kp=''!
else
        echo "arc compute undetermined"
        exit=0
fi



# clean up the pdb in /mnt/c/Users/james/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe
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
echo "remove hydrogens
remove not alt ''+A
alter all, alt=''
save $name.clean.pdb, $name, state=0" >> clean.pml
/mnt/c/Users/james/AppData/Local/Schrodinger/PyMOL2/PyMOLWin.exe -qc clean.pml
sleep 5
grep -v ANISOU $name.clean.pdb > $name.cleaner.pdb

sed -i 's/HIE/HIS/g;s/HSD/HIS/g;s/HSE/HIS/g;s/CYX/CYS/g' $name.cleaner.pdb
sed -i 's/CD  ILE/CD1 ILE/g;s/OT1/O  /g;s/OT2/OXT/g' $name.cleaner.pdb

# if its this big then we really ought to use high memory
highmem=no
if [ $(wc -l < $name.cleaner.pdb) -gt 500000 ] ; then
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
                        sed -i 's/./0/73;s/./0/74;s/./0/75;s/./ /76' $name.cleaner.pdb
                fi

        fi
        for j in $(sed -n "/^$MODEL/,/^ENDMDL/{p;/^ENDMDL/q}" $name.cleaner.pdb | cut -c 73-76 | grep "\S" | sort | uniq ) ; do
                segi=$(echo 00$count | rev | cut -c 1-4 | rev )
                sed -n "/^$MODEL/,/^ENDMDL/{p;/^ENDMDL/q}" $name.cleaner.pdb | grep "     $j " | sed "s/      $j /      $segi /g" >> $name.segment.pdb
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
awk '/[A-Z][A-Z][A-Z ] '$j' .* '$i' / && !x {print "TER"; x=1} 1' $name.TER.pdb > temp.pdb
mv temp.pdb $name.TER.pdb
fi
done
done
echo 'TER' >> $name.TER.pdb
uniq $name.TER.pdb > $name.TER2.pdb

# check for glycam residues and add TER
grep ATOM $name.TER2.pdb | cut -c 18-20 | sort | uniq > resn-list.tx
cat /home/james/anaconda3/pkgs/ambertools-19-h0d7ec52_0/dat/leap/prep/GLYCAM_06j-1.prep | grep INT | cut -c 1-3 > glycam-list.tx

for sugar in $(grep -f resn-list.tx glycam-list.tx ) ; do 
grep $sugar $name.TER2.pdb | cut -c 18-26 | sort | uniq > glycam2-list.tx
for glyresn in $(seq 1 $(wc -l < glycam2-list.tx) ) ; do 
glyres=$(sed -n ''$glyresn'p' glycam2-list.tx)
awk '/'"$(grep "$glyres" $name.TER2.pdb | head -1 | cut -c 1-26)"'/ && !x {print "TER"; x=1} 1 ' $name.TER2.pdb > temp.pdb
echo adding $sugar $glyres 
#awk '/'$glyres'/ && !x {print "TER"; x=1} 1' $name.TER.pdb > temp.pdb
mv temp.pdb $name.TER2.pdb
done
done



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

# add disulfide bond fix. here we change residues to CYX and add conect records.

for i in $( seq 1 $( wc -l < ssbond.tx ) ) ; do 

sed -i "/$(sed -n ''$i'p'  ssbond.tx | grep -o  "CYS .....[0-9]" | head -1 | cut -c 5,7-)/s/CYS/CYX/g" $name.ssbond.pdb
sed -i "/$(sed -n ''$i'p'  ssbond.tx | grep -o  "CYS .....[0-9]" | tail -1 | cut -c 5,7-)/s/CYS/CYX/g" $name.ssbond.pdb

echo "CONECT "$(grep "$(sed -n ''$i'p'  ssbond.tx | grep -o  "CYS .....[0-9] " | head -1 | cut -c 5,7- )" $name.ssbond.pdb | grep SG | awk '{print $2}')" "$(grep "$(sed -n ''$i'p'  ssbond.tx | grep -o  "CYS .....[0-9]" | tail -1 | cut -c 5,7-)" $name.ssbond.pdb | grep SG | awk '{print $2}') >> $name.ssbond.pdb
done

# pass to amber to check through structure.
# pdb4amber is broken for more than 9999 residues
if [ $( grep ' CA ' $name.ssbond.pdb | wc -l ) -lt 9999 ] ; then
        #pdb4amber -p $name.ssbond.pdb > $name.amber.pdb
        cp $name.ssbond.pdb $name.amber.pdb
else
        cp $name.ssbond.pdb $name.amber.pdb
fi


# failure-cleanup
# if the run fails for any reason, run this script to return to a native folder (does not remove pdbs)
echo "
rm min* eqi* pro*
rm -r $outdir
rm *.tx *clean.pdb *cleaner.pdb *segment.pdb *TER.pdb *TER2.pdb *ssbond.pdb 
rm clean.pml dmp.py distance-matrix* 
" > failure-cleanup.sh
chmod 755 failure-cleanup.sh


########################################################################################
# BODY TEXT - EXPLICIT SOLVATON
########################################################################################

# set boundary for edge waters
gap=14.0

# Generate TEST tleap
rm leap1.log 2>/dev/null
echo "source leaprc.GLYCAM_06j-1
source leaprc.protein.ff14SB
source leaprc.water.tip3p

mol = loadPdb $name.amber.pdb
solvateoct mol TIP3PBOX $gap iso
charge mol
quit" > tleap.test.in

tleap -f tleap.test.in > leap1.log

# get waters to calculate ions
waters=$(grep Added leap1.log | tail -1 | awk '{print $2}')

# get volume from min/max dimentions for salt concentration.
calcvolume=$(grep Volume leap1.log | awk '{print $2}' | awk -F'.' '{print $1}')

# adjust volume by subtracting protein atoms - arbitratry factor of 5!
advolume=$(( $(tail $name.amber.pdb | awk '{print $2}' | grep "\S" | tail -1) * 10 ))
volume=$(( $calcvolume - $advolume ))

# find 150 mM salt ions
ionsold=$(( $volume * 903 / 10000000 ))

# find number of ions from water
ions=$(echo "(($waters * 0.0187 * $saltconc) + 0.5) / 1" | bc)

echo "ions old = $ionold, ions new = $ions"

# find charge to neutrilise
charge=$(grep "Total unperturbed charge" leap1.log | awk '{print $4}' | awk -F'.' '{print $1}')
if [ $charge -gt 0 ] ; then

        nions=$(( $charge + $ions ))
        pions=$ions
elif [ $charge -lt 0 ] ; then

        pions=$(( $(echo $charge | cut -c 2- ) + $ions ))
        #pions=$(( $ncharge + $ions ))
        nions=$ions
else
        pions=$ions
        nions=$ions
fi

# Generate solvated volume
rm leap.log 2>/dev/null
echo "source leaprc.GLYCAM_06j-1
source leaprc.protein.ff14SB
source leaprc.water.tip3p

mol = loadPdb $name.amber.pdb
solvateoct mol TIP3PBOX $gap iso
addIonsRand mol Na+ $pions Cl- $nions
savepdb mol $name.explicit.pdb
saveAmberParm mol $name.explicit.parm7 $name.explicit.rst7
quit" > tleap.solvate.in

# With the ff19SB force field there is an updated water model.
# can use : source leaprc.protein.ff19SB
# with    : source leaprc.water.opc
#                 : loadamberparams frcmod.ions1lm_126_hfe_opc
#                 : solvateOct ramp SPCBOX 14.0

tleap -f tleap.solvate.in

avolume=$(grep Volume leap.log | awk '{print $2}' | awk -F'.' '{print $1}')
mv leap.log leap.solvate.log
rm tleap.in 2>/dev/null

# compare calclutated volumes
echo "EstiVolume = $volume
CalcVolume = $avolume
Vol-Change = $(( $volume - $avolume ))"

echo "hmassrepartition
outparm $name.explicit.hmr.parm7
quit" > parm.in
parmed -i parm.in $name.explicit.parm7

# SIMULATION ARC JOBS (DEFAULTS)
########################################################################################

# set arc parameters
warc=$(echo $arc | cut -c 1)
if [ $warc = 3 ] ; then
        COMPUTE=arc3.leeds.ac.uk
        CPU=24
        kp=$(echo $arc | cut -c 2)
        if [ $kp = p ] ; then
                GPU=p100
        else
                GPU=k80
        fi

elif [ $warc = 4 ] ; then
        COMPUTE=arc4.leeds.ac.uk
        CPU=40
        GPU=v100
else
        echo "arc compute undetermined"
        exit=0ACO_AC.lib
fi

# list jobs
ls -v minq* > joblist.tx
ls -v eqiq* >> joblist.tx
ls -v proq* >> joblist.tx

first=1
for i in $(cat joblist.tx) ; do
        nsub=$(echo $i | sed 's/.mdin//g')
        if [ $(echo $nsub | cut -c 1-3 ) = min ] ; then
                echo "#!/bin/bash
#$ -pe ib $CPU
#$ -l h_vmem=4.5G
#$ -l h_rt=48:0:0
#$ -cwd -V
#$ -m abe
#$ -M chmjro@leeds.ac.uk
#$ -N $nsub-$rand
#$ -hold_jid $last
module add amber
./$nsub" > $nsub-$rand.sh
                if [ $first = 1 ] ; then
                        sed -i '/old_jid/d' $nsub-$rand.sh
                        first=2
                fi

        else
                echo "#!/bin/bash
#$ -l coproc_$GPU=1
#$ -l h_rt=48:0:0
#$ -cwd -V
#$ -m abe
#$ -M chmjro@leeds.ac.uk
#$ -N $nsub-$rand
#$ -hold_jid $last
module add amber/20gpu
./$nsub" > $nsub-$rand.sh
        fi
        last=$nsub-$rand
done

# check for highmem requirement
if [ $highmem = yes ] ; then
        ls -v min*$rand.sh > joblist2.tx
        for i in $(cat joblist2.tx) ; do ACO_AC.lib
                sed -i "/^#!/bin/bash/a  node_type=$CPU\core-768G" $i
                sed -i "s/vmem=4.5G/vmem=18G/g" $i
        done
fi

#

        # . . . with md input files
cp min*.mdin $outdir/.
cp eqi*.mdin $outdir/.
cp pro*.mdin $outdir/.

        # . . . with structure files
mv $name.explicit.hmr.parm7 $outdir/initial.parm7
mv $name.explicit.rst7 $outdir/initial.rst7
mv $name.explicit.pdb $outdir/initial.pdb
rm *.explicit.parm7

        # . . . with md submission files
cp *$rand.sh $outdir/.
cp minq* $outdir/.
cp eqiq* $outdir/.
cp proq* $outdir/.

chmod 755 $outdir/*

        # pass to arc and submit jobs
ssh $PROXY $USER@$COMPUTE mkdir $remotedir/$date-amber-$name
scp $PROXY -r $date-amber-$name $USER@$COMPUTE:$remotedir/. 1>/dev/null
for i in $(ls -v min*$rand.sh) ; do
ssh $PROXY $USER@$COMPUTE "cd $remotedir/$date-amber-$name ; qsub $i" > $date.tx
done
for i in $(ls -v eqi*$rand.sh) ; do
ssh $PROXY $USER@$COMPUTE "cd $remotedir/$date-amber-$name ; qsub $i" > $date.tx
done
for i in $(ls -v pro*$rand.sh) ; do
ssh $PROXY $USER@$COMPUTE "cd $remotedir/$date-amber-$name ; qsub $i" > $date.tx
done

rm min* eqi* pro*


# POST Production
########################################################################################
# run this script once complete!!

echo "#!/bin/bash

ls -v pro*.nc > nc.list

# determine the fixed positions for auto imaging
epro=\$(grep ' CA  [A-Z][A-Z][A-Z] '  initial.pdb | nl | tail -1 | awk '{print \$1}')

echo 'parm initial.parm7
trajin eqi1.nc' > cpptraj.in
sed 's/^/trajin /g' nc.list > nct.list
cat nct.list >> cpptraj.in

echo 'strip :WAT outprefix traj
center :1-'\$epro' mass origin
unwrap :1-'\$epro' center
trajout traj.dcd
trajout traj.pdb onlyframes 1
trajout traj.rst7 onlyframes 1
cluster clusters 5 repout topcluster repfmt pdb out cnumvtime.dat summary avg.summary.dat
readdata pro*.mdout name MDOUT
writedata temp.dat MDOUT[TEMP] time 0.00002
writedata etot.dat MDOUT[Etot] time 0.00002
writedata density.dat MDOUT[Density] time 0.00002
writedata volume.dat MDOUT[VOLUME] time 0.00002
2drms :1-'\$epro'@CA&!@H= rmsout 2drms-all-residues.gnu
average average.pdb pdb
rms ToproFirst :1-'\$epro'&!@H= first out rmsdpro.agr mass
atomicfluct :1-'\$epro'&!@H= byres out fluctpro.agr
atomicfluct :1-'\$epro'@CA&!@H= byres bfactor out fluctbfact.agr
run
exit ' >> cpptraj.in

module add amber/16
cpptraj -i cpptraj.in

mkdir trajout
mv traj* trajout/.
cp *.out trajout/.
cp *info trajout/.
mv *.dat trajout/.
mv topcluster* trajout/.
mv *agr trajout/.
mv average.pdb trajout/.
mv avg.summary.dat trajout/.
sed -i \"s/pause -1/set terminal png size 1200,1200\nset output 'rmsd-matrix.png'\nreplot\nquit/g\" 2drms-all-residues.gnu
mv 2drms-all-residues.gnu trajout/.
dmp.py
" > jr-amber-post-process.sh
chmod 755 jr-amber-post-process.sh
scp $PROXY -r jr-amber-post-process.sh $USER@$COMPUTE:$remotedir/$date-amber-$name/. 1>/dev/null

#cleanup
rm *.tx *clean.pdb *cleaner.pdb *segment.pdb *TER.pdb *TER2.pdb *ssbond.pdb *amber.pdb 
rm clean.pml dmp.py distance-matrix* 

exit 0
