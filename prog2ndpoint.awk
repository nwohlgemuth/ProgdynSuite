BEGIN {
# aug 2013 includes molecular rotation, ability to do multiple NMR calculations, ONIOM with link atoms, 
#nonstandard routes, checks more kinds of energies at point 2
#Aug 2010 added etolerance to make it controllable from progdyn.conf, made it so that DRP does not check energy
# aug 2008 added to atom list so handles 1 to 17 without change needed
# version Feb 2008 incorporates methodfile, boxon and boxsize, though this point unaffected by box
# version Jan 2008 - allows for ONIOM jobs, fixed atoms
# version Sept 9, 2005 - incorportates meth3, meth4, meth5, meth6, but not yet rotation
# read progdyn.conf for configuration info

# default parameters, including quassiclassical, no displacements, transition state, not a DRP
# do not change these - rather, change progdyn.conf to set the parameters
initialDis=0; timestep=1E-15; scaling=1.0; temp=298.15
classical=0; numimag=1; DRP=0; cannonball=0
memory=20000000
diag=1; checkpoint="g09.chk"; searchdir="positive"; boxon=0
boxsize=10; maxAtomMove=0.1; title1="you"; title2="need"
title3="a"; title4="progdyn.conf"; processors=1; highlevel=999; linkatoms=0
etolerance=1
geometry="nonlinear";nonstandard=0
NMRtype=none;NMRevery=9999999

#initialization
i=1;j=1;k=1
c=29979245800; h=6.626075E-34; avNum=6.0221415E23
RgasK=0.00198588; RgasJ=8.31447
numAtoms=0; atomnumber=0
getline < "runpointnumber"
runpointnum = $1

blankLineTester=10
while (blankLineTester>1) {
   getline < "progdyn.conf"
   if ($1=="method") method=$2
   if ($1=="method2") meth2=$2
   if ($1=="charge") charge=$2
   if ($1=="multiplicity") multiplicity=$2
   if ($1=="memory") memory=$2
   if ($1=="processors") processors=$2
   if ($1=="checkpoint") checkpoint=$2
   if ($1=="timestep") timestep=$2
   if ($1=="diagnostics") diag=$2
   if ($1=="method3") meth3=$2
   if ($1=="method4") meth4=$2
   if ($1=="method5") meth5=$2
   if ($1=="method6") meth6=$2
   if ($1=="method7") meth7=$2
   if ($1=="highlevel") highlevel=$2
   if ($1=="linkatoms") linkatoms=$2
   if ($1=="fixedatom1") fixedatom1=$2
   if ($1=="fixedatom2") fixedatom2=$2
   if ($1=="fixedatom3") fixedatom3=$2
   if ($1=="fixedatom4") fixedatom4=$2
   if ($1=="DRP") DRP=$2
   if ($1=="methodfile") methodfilelines=$2
   if ($1=="killcheck") killcheck=$2
   if ($1=="etolerance") etolerance=$2
   if ($1=="reversetraj") reversetraj=$2
   if ($1=="NMRmethod") nmrmethod=$2
   if ($1=="NMRmethod2") nmrmethod2=$2
   if ($1=="NMRmethod3") nmrmethod3=$2
   if ($1=="NMRtype") nmrtype=$2
   if ($1=="NMRevery") nmrevery=$2
   if ($1=="nonstandard") nonstandard=$2
   if ($1=="title") {
      title1=$2
      title2=$3
      title3=$4
      title4=$5
      }
   blankLineTester=length($0)
   }

if (diag>=1) print "***************** starting prog2ndpoint *****************" >> "diagnostics"
if (diag>=1) print "method,charge,multiplicity,memory" >> "diagnostics"
if (diag>=1) print method,charge,multiplicity,memory >> "diagnostics"
if (diag>=1) print "processors,checkpoint,title" >> "diagnostics"
if (diag>=1) print processors,checkpoint,title1,title2,title3,title4 >> "diagnostics"

#get the isomer number from file
getline < "isomernumber"
isomernum = $1

#get forward or reverse from skipstart if it exists
getline < "skipstart"
trajdirection = $1

print "%nproc=" processors
print "%mem=" memory
if (killcheck!=1) print "%chk=" checkpoint
if (nonstandard==0) {
   print "# " method " force scf=(tight,nosym) "
   if (meth2=="unrestricted") print "guess=mix" #for unrestricted calculations
   if (meth2=="read") print "guess=tcheck" #for reading orbitals from check, sometimes faster, sometimes not
   if (length(meth3)>2) print meth3
   if (length(meth4)>2) print meth4
   }
if (nonstandard==1) {
   print "# "
   print "nonstd"
   system("cat nonstandard")
   }
print ""
print title1,title2,title3,title4
print "runpoint ",runpointnum
print "runisomer ", isomernum
print ""
print charge,multiplicity

# ok, now we have to figure the second point.  this should be 
# x(t) = x + v*t + 1/2*F*t^2/m
# so we need to set up arrays for position, velocity, and force

#read in number of atoms, geometry, masses, and velocity from geoPlusVel
getline < "geoPlusVel"
numAtoms=$1
# geometry
for (i=1;i<=numAtoms;i++) {
   getline < "geoPlusVel" 
   weight[i]=$5
   atSym[i]=$1
   for (j=1;j<=3;j++) {
      geoArr[i,j]=$(1+j)
      }
   }
#velocities
for (i=1;i<=numAtoms;i++) {
   getline < "geoPlusVel"
   for (j=1;j<=3;j++) {
      velArr[i,j]=$j
      }
   }

#now we go ahead and add the velocities
for (i=1;i<=numAtoms;i++) {
   for (j=1;j<=3;j++) {
      arr[i,j]=velArr[i,j]+geoArr[i,j]
      if (trajdirection=="reverserestart") arr[i,j]=geoArr[i,j]-velArr[i,j]
      }
   if ((diag>1) && (i==1)) print "geometry after adding velocities" >> "diagnostics"
   if (diag>1) print arr[i,1],arr[i,2],arr[i,3] >> "diagnostics"
   }

#pull out other information useful for testing whether total energy is right or bad
blankLineTester=10
while (blankLineTester>1) {
   getline < "geoPlusVel"
   if ($4=="desired=") desiredModeEnK=$5
   if ($4=="modes=") {
      KEinitmodes=$5
      KEinittotal=$9
      }
   if ($11=="potential") potentialE=$13
   blankLineTester=length($0)
   }
#get initial geometry into file traj
print numAtoms >> "traj"
print potentialE,title1,title2,title3,title4,"runpoint 1 ","runisomer ",isomernum >> "traj"
for (i=1;i<=numAtoms;i++) {
   print atSym[i],geoArr[i,1],geoArr[i,2],geoArr[i,3] >> "traj"
   }
#added by Samae on 102910
scfcount=0
} # end of BEGIN

#pull out the potential energy
/SCF Done/ || /EUMP2 =/ || / Energy=/ {
if (($1=="Energy=") && ($3=="NIter="))  newPotentialE=$2
if (($1=="SCF") && (scfcount==0)) newPotentialE=$5
if ($1=="E2") {
   tempstring=$6
   split(tempstring, arr10, "D")
   newPotentialE=arr10[1]*(10^arr10[2])
   }
newPotentialEK=(newPotentialE-potentialE)*627.509
if ($1=="SCF") {
   if (scfcount==0) {
      pddga=$5
      }
   if (scfcount==1) {
      qm=$5
      }
   if (scfcount==2) {
      pddgb=$5
      pddgc=(pddga-pddgb)
      newPotentialE=(qm+pddgc)
      newPotentialEK=(newPotentialE-potentialE)*627.509
      }
   scfcount++
   }
}

# now we go ahead and translate the forces and add them
(/        1    / || /        2    / || /        3    / || /        4    / || /        5    / || /        6    / || /        7    / || /        8    / || /        9    / || /       10    / || /       11    / || /       12    / || /       13    / || /       14    / || /       15    / || /       16    / || /       17    / || /       18    / || /       19    / || /       20    / || /       21    / || /       22    / || /       23    / || /       24    / || /       25    / || /       26    / || /       27    / || /       28    / || /       29    / || /       30    / || /       31    / || /       32    / || /       33    / || /       34    / || /       35    /) && length($3) > 9 {
i=$1
for (j=1;j<=3;j++) {
   forceArr[i,j]=$(2+j)    #the raw units of the forces are Hartree/Bohr
   }
if ((diag>1) && (i==1)) print "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" >> "diagnostics"
if (diag>1) print i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3] >> "diagnostics"
}

END {
#put out Echeck but only if not a DRP
if (DRP==0) {
   print "trajectory #",isomernum >> "Echeck"
   print "point 1 potential E=",newPotentialEK,"   point 1 kinetic E=",KEinitmodes,"  Total=",newPotentialEK+KEinitmodes >> "Echeck"
   print "desired total energy=", desiredModeEnK >> "Echeck"
   if ((newPotentialEK+KEinitmodes)>(desiredModeEnK+etolerance)) print "XXXX bad total Energy" >> "Echeck"
   if ((newPotentialEK+KEinitmodes)<(desiredModeEnK-etolerance)) print "XXXX bad total Energy" >> "Echeck"
   }
# turn the forces into motion
for (i=1;i<=numAtoms;i++) {
   for (j=1;j<=3;j++) {
# conversions here take force to J/angstrom, 1E20 converts to kg angstroms / s^2, then mult time (s^s) and divide by weight in kg to get angstroms
      forceArr[i,j]=0.5*1E20*forceArr[i,j]*627.509*(4184/(0.529177*avNum))*(timestep^2)/(weight[i]/(avNum*1000))
# for simplicity, DRPs will throw away the forces at the second pont.  This means that if we are not at a saddlepoint, point 2 = point 1 but this is a minor waste
      if (DRP==1) forceArr[i,j]=0
      arr[i,j]=arr[i,j]+forceArr[i,j]
# if atoms are fixed, replace calcd new position by original position
      if ((i==fixedatom1) || (i==fixedatom2) || (i==fixedatom3) || (i==fixedatom4)) arr[i,j]=geoArr[i,j]
      }
   if ((diag>1) && (i==1)) print "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" >> "diagnostics"
   if (diag>1) print i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3] >> "diagnostics"
   printf("%s %.7f %.7f %.7f",atSym[i],arr[i,1],arr[i,2],arr[i,3])
   if ((i>highlevel) && (i<=highlevel+linkatoms)) printf(" %s","M H")
   if (i>(highlevel+linkatoms)) printf(" %s","M")
   print ""
   }
print ""
if (length(meth5)>2) print meth5
if (length(meth6)>2) print meth6
if (methodfilelines>=1) {
   for (i=1;i<=methodfilelines;i++) { 
      getline < "methodfile" 
      print $0
      }
   }
if ((nmrtype>0) && ((runpointnum % nmrevery)==0)) {
   print "--link1--"
   print "%nproc=" processors
   print "%mem=" memory
   print "%chk=" checkpoint
   print "# " nmrmethod " nmr=giao geom=check"
   if (nmrmethod==method) print "guess=tcheck"
   if (length(meth7)>2) print meth7
   print ""
   print title1,title2,title3,title4
   print "runpoint ",runpointnum
   print "runisomer ",isomernum
   print ""
   print charge,multiplicity
   }
print ""
if ((nmrtype>1) && ((runpointnum % nmrevery)==0)) {
   print "--link1--"
   print "%nproc=" processors
   print "%mem=" memory
   print "%chk=" checkpoint
   print "# " nmrmethod2 " nmr=giao geom=check"
   if (length(meth7)>2) print meth7
   print ""
   print title1,title2,title3,title4
   print "runpoint ",runpointnum
   print "runisomer ",isomernum
   print ""
   print charge,multiplicity
   }
print ""
if ((nmrtype>2) && ((runpointnum % nmrevery)==0)) {
   print "--link1--"
   print "%nproc=" processors
   print "%mem=" memory
   print "%chk=" checkpoint
   print "# " nmrmethod3 " nmr=giao geom=check"
   if (length(meth7)>2) print meth7
   print ""
   print title1,title2,title3,title4
   print "runpoint ",runpointnum
   print "runisomer ",isomernum
   print ""
   print charge,multiplicity
   }
print ""
#get second geometry into file traj
print numAtoms >> "traj"
print newPotentialE,title1,title2,title3,title4,"runpoint ",runpointnum,"runisomer ",isomernum >> "traj"
for (i=1;i<=numAtoms;i++) {
   print atSym[i],arr[i,1],arr[i,2],arr[i,3] >> "traj"
   }
}

