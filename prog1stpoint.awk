BEGIN {
# aug 2013 includes molecular rotation, ability to do multiple NMR calculations, ONIOM with link atoms, 
#nonstandard routes, 
# aug 2010 changed so that it is more careful in reading in from geoPlusVel
#     removed some default parameters that should always be defined
# Jan 2009 - a number of little changes to improve reporting, precision, etc
# aug 2008 added to atom list so handles H to Cl without change needed
# version Feb 2008 incorporates methodfile, boxon and boxsize, though this point unaffected by box
# version Jan 2008 - allows for ONIOM jobs, fixed atoms
# version Sept 2005 - incorportates meth3, meth4, meth5, meth6, but not yet rotation
# this program creates the first input file for g09
# the title should be changed as appropriate
# the isomer number comes from a file isomernumber

# default parameters, including quassiclassical, no displacements, transition state, not a DRP
# do not change these - rather, change progdyn.conf to set the parameters
initialDis=0; timestep=1E-15; scaling=1.0; temp=298.15
classical=0; numimag=1; DRP=0; cannonball=0
memory=20000000
diag=1; checkpoint="g09.chk"; searchdir="positive"; boxon=0
boxsize=10; maxAtomMove=0.1; title1="you"; title2="need"
title3="a"; title4="progdyn.conf"; processors=1; highlevel=999; linkatoms=0
geometry="nonlinear";nonstandard=0
nmrtype=0;nmrevery=9999999

#initialization
i=1;j=1;k=1
c=29979245800; h=6.626075E-34; avNum=6.0221415E23
RgasK=0.00198588; RgasJ=8.31447
numAtoms=0; atomnumber=0
getline < "runpointnumber"
runpointnum = $1

# read progdyn.conf for configuration info
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
   if ($1=="methodfile") methodfilelines=$2
   if ($1=="killcheck") killcheck=$2
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

if (diag==1) print "***************** starting prog1stpoint *****************" >> "diagnostics"
if (diag==1) print "method,charge,multiplicity,memory" >> "diagnostics"
if (diag==1) print method,charge,multiplicity,memory >> "diagnostics"
if (diag==1) print "processors,checkpoint,title" >> "diagnostics"
if (diag==1) print processors,checkpoint,title1,title2,title3,title4 >> "diagnostics"

getline < "isomernumber"
isomernum = $1
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
#velocities not needed for 1st point
for (i=1;i<=numAtoms;i++) {
   getline < "geoPlusVel"
   for (j=1;j<=3;j++) {
      velArr[i,j]=$j
      }
   }

print "%nproc=" processors
print "%mem=" memory
if (killcheck!=1) print "%chk=" checkpoint
if (nonstandard==0) {
   print "# " method " force scf=(tight,nosym) "
   if (meth2=="unrestricted") print "guess=mix" #for unrestricted calculations
   if (length(meth3)>2) print meth3
   if (length(meth4)>2) print meth4
   }
if (nonstandard==1) {
   print "# "
   print "nonstd"
   system("cat nonstandard")
   }
print ""
# make the title four words exactly, leaving out spaces if necessary
print title1,title2,title3,title4
print "runpoint ",runpointnum
print "runisomer ", isomernum
print ""
print charge,multiplicity
}

END {
for (i=1;i<=numAtoms;i++) {
   printf("%s %.7f %.7f %.7f",atSym[i],geoArr[i,1],geoArr[i,2],geoArr[i,3])
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
}

