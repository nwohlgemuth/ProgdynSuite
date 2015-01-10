BEGIN { #this is the main routine for generating new .com files by the Verlet algorithym 
# aug 2013 includes molecular rotation, ability to do multiple NMR calculations, ONIOM with link atoms, 
#nonstandard routes, monitoring of loads and randomized NMR calcs, CCSD(T) nmr calculations, making a ZMAT file for CFOUR
# Aug 2010 increased elements handled automatically but only up to bromine!
# Jan 2009 - a number of little changes to improve reporting, precision, etc
# Nov 2008 added ability to handle DRPs
# Aug 2008 added long list of atoms to handle 1-17 without change
# May 2008 added option to put out velocities in vellist - make diag=3
# version Feb 2008 incorporates methodfile, boxon and boxsize
# version Jan 2008 incorporates fixed atoms, oniom, and velocity damping
# version August 2007 incorporates keepevery to decrease size of dyn file
# version Sept 11, 2005 - incorportates meth3, meth4, meth5, meth6, but not yet rotation

# default parameters, including quassiclassical, no displacements, transition state, not a DRP
# do not change these - rather, change progdyn.conf to set the parameters
initialDis=0; timestep=1E-15; scaling=1.0; temp=298.15
classical=0; numimag=1; DRP=0; cannonball=0
memory=20000000
diag=1; checkpoint="g09.chk"; searchdir="positive"; boxon=0
boxsize=10; maxAtomMove=0.1; title1="you"; title2="need"
title3="a"; title4="progdyn.conf"; processors=1; highlevel=999; linkatoms=0
damping=1;nonstandard=0
nmrtype=0;nmrevery=9999999;nmrcc=0;nmrrand=0;nmrdo=0
thermostat=0;thermostatemult=1.00

#initialization
srand(PROCINFO["pid"])
i=1;j=1;k=1
c=29979245800; h=6.626075E-34; avNum=6.0221415E23
RgasK=0.00198588; RgasJ=8.31447
numAtoms=0; atomnumber=0
conver1=4.184E26 #dividing by this converts amu angs^2 /s^2 to kcal/mol
OFS="     "

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
   if ($1=="temperature") temp=$2
   if ($1=="thermostat") thermostat=$2
   if ($1=="thermostatmult") thermostatmult=$2
   if (thermostatmult>1) thermostatmult=1/thermostatmult
   if ($1=="method3") meth3=$2
   if ($1=="method4") meth4=$2
   if ($1=="method5") meth5=$2
   if ($1=="method6") meth6=$2
   if ($1=="method7") meth7=$2
   if ($1=="highlevel") highlevel=$2
   if ($1=="linkatoms") linkatoms=$2
   if ($1=="keepevery") keepevery=$2
   if ($1=="fixedatom1") fixedatom1=$2
   if ($1=="fixedatom2") fixedatom2=$2
   if ($1=="fixedatom3") fixedatom3=$2
   if ($1=="fixedatom4") fixedatom4=$2
   if ($1=="boxon") boxon=$2
   if ($1=="boxsize") boxsize=$2
   if ($1=="DRP") DRP=$2
   if ($1=="maxAtomMove") maxAtomMove=$2
   if ($1=="methodfile") methodfilelines=$2
   if ($1=="killcheck") killcheck=$2
   if ($1=="damping") damping=$2
   if ($1=="NMRmethod") nmrmethod=$2
   if ($1=="NMRmethod2") nmrmethod2=$2
   if ($1=="NMRmethod3") nmrmethod3=$2
   if ($1=="NMRtype") nmrtype=$2
   if ($1=="NMRevery") nmrevery=$2
   if ($1=="NMRrand") nmrrand=$2
   if ($1=="loadlimit") loadlimit=$2
   if ($1=="NMRcc") nmrcc=$2
   if ($1=="nonstandard") nonstandard=$2
   if ($1=="title") {
      title1=$2
      title2=$3
      title3=$4
      title4=$5
      }
   blankLineTester=length($0)
   }

if (diag>=1) print "***************** starting progdynb *****************" >> "diagnostics"
if (diag>=1) print "method,charge,multiplicity,memory" >> "diagnostics"
if (diag>=1) print method,charge,multiplicity,memory >> "diagnostics"
if (diag>=1) print "processors,checkpoint,title" >> "diagnostics"
if (diag>=1) print processors,checkpoint,title1,title2,title3,title4 >> "diagnostics"

# get number of atoms and weights from geoPlusVel, and previous geometries from old and older
getline < "geoPlusVel"
numAtoms=$1
for (i=1;i<=numAtoms;i++) {
   getline < "geoPlusVel" 
   weight[i]=$5; atSym[i]=$1
   }

for (at=1;at<=numAtoms;at++) {
   getline < "old"
   oldarr[at,1]=$4; oldarr[at,2]=$5; oldarr[at,3]=$6
   }

for (at=1;at<=numAtoms;at++) {   
   getline < "older"   
   olderarr[at,1]=$4; olderarr[at,2]=$5; olderarr[at,3]=$6
   }   

#for DRPs read in oldAdjForces and maxAtomMove
if (DRP==1) {
   for (at=1;at<=numAtoms;at++) {
      getline < "oldAdjForces"
      oldForce[at,1]=$1; oldForce[at,2]=$2; oldForce[at,3]=$3
      }
   getline < "maxMove"
   if (($1<maxAtomMove) && ($1>0)) maxAtomMove=$1
   if (maxAtomMove<0.000001) maxAtomMove=0.000001
   }

# record atom velocities for IVR analysis.  This is actually the velocity in the previous run, which is the easiest to calculate.
getline < "isomernumber"
isomernum = $1
getline < "runpointnumber"
runpointnum = $1
# routine to control whether NMR calculations are done.  
if ((nmrrand==0) && ((runpointnum % nmrevery)==0)) nmrdo=1
if ((nmrrand==1) && (rand()<(1/nmrevery))) nmrdo=1
getline < "uptimelist"
x=1.0001*substr($10,1,3);if (x<8) x=8
# turn of nmrs if load is too high - this is under control of loadlimit parameter in progdyn.conf and requires proganal to make uptimelist
if ((nmrrand==1) && (x>loadlimit)) nmrdo=0

if (diag==3) print "runpoint ",runpointnum-1,"runisomer ",isomernum >> "vellist"
for (at=1;at<=numAtoms;at++) {
   atomVel=((oldarr[at,1]-olderarr[at,1])^2 + (oldarr[at,2]-olderarr[at,2])^2 +(oldarr[at,3]-olderarr[at,3])^2)^.5
   KEatomstotal=KEatomstotal+0.5*weight[at]*(atomVel^2)/((timestep^2)*conver1)
   if (diag==3) print atomVel >> "vellist"
   }
apparentTemp=KEatomstotal*2/(3*RgasK*numAtoms)
if (diag==4) print "KEatomstotal",KEatomstotal,"apparent Temperature",apparentTemp >> "vellist"
if (thermostat==1) {
   if (diag<4) print "KEatomstotal",KEatomstotal,"desired temperature",temp,"apparent Temperature",apparentTemp >> "vellist"
   if (apparentTemp>temp) damping=thermostatmult
   if (apparentTemp<temp) damping=1/thermostatmult
   }
}

#pull out the potential energy
/SCF Done/ || /EUMP2 =/ || / Energy=/ || /ONIOM:/ {
if (($1=="Energy=") && ($3=="NIter="))  newPotentialE=$2
if ($1=="SCF") newPotentialE=$5
if ($2=="extrapolated")  newPotentialE=$5
if ($1=="E2") {
   tempstring=$6
   split(tempstring, arr10, "D")
   newPotentialE=arr10[1]*(10^arr10[2])
   }
}

#must adjust next line for weird atoms
(/        1    / || /        2    / || /        3    / || /        4    / || /        5    / || /        6    / || /        7    / || /        8    / || /        9    / || /       10    / || /       11    / || /       12    / || /       13    / || /       14    / || /       15    / || /       16    / || /       17    / || /       18    / || /       19    / || /       20    / || /       21    / || /       22    / || /       23    / || /       24    / || /       25    / || /       26    / || /       27    / || /       28    / || /       29    / || /       30    / || /       31    / || /       32    / || /       33    / || /       34    / || /       35    /) && length($3) > 9 {
i=$1
for (j=1;j<=3;j++) {
   forceArr[i,j]=$(2+j)    #the raw units of the forces are Hartree/Bohr
   }
if ((diag>1) && (i==1)) print "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" >> "diagnostics"
if (diag>1) print i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3] >> "diagnostics"
}

END {
#############routine for DRPs##############
if (DRP==1) {
   maxForce=0;oscillTest=0
   for (i=1;i<=numAtoms;i++) {
      for (j=1;j<=3;j++) {
# conversions here take force to J/angstrom, 1E20 converts to kg angstroms / s^2, then mult time (s^s) and divide by weight in kg to get angstroms
         forceArr[i,j]=1E20*forceArr[i,j]*627.509*(4184/(0.529177*avNum))*(timestep^2)/(weight[i]/(avNum*1000))
         oscillTest=oscillTest+forceArr[i,j]*oldForce[i,j]
         if (forceArr[i,j]>maxForce) maxForce=forceArr[i,j]
         if ((0-forceArr[i,j])>maxForce) maxForce=-forceArr[i,j]
         }
      if (i==1) printf("% .8f % .8f % .8f \n",forceArr[1,1],forceArr[1,2],forceArr[1,3])  > "oldAdjForces"
      if (i>1) printf("% .8f % .8f % .8f \n",forceArr[i,1],forceArr[i,2],forceArr[i,3])  >> "oldAdjForces"
      }
   print "oscillTest ",oscillTest >> "oldAdjForces"
   if (oscillTest<0) {
      maxAtomMove = maxAtomMove*0.5
      print maxAtomMove > "maxMove"
      }
   if (oscillTest>0) {
      maxAtomMove = maxAtomMove*1.2
      print maxAtomMove > "maxMove"
      }
   print "maxAtomMove ",maxAtomMove >> "oldAdjForces"
   forceMult=maxAtomMove/maxForce
   for (i=1;i<=numAtoms;i++) {
      for (j=1;j<=3;j++) {
         newarr[i,j]=oldarr[i,j]+forceMult*forceArr[i,j]
         }
      }
   }
########

#############normal routine for Verlet ##############
if (DRP==0) {
   for (i=1;i<=numAtoms;i++) {
      for (j=1;j<=3;j++) {
# conversions here take force to J/angstrom, 1E20 converts to kg angstroms / s^2, then mult time (s^s) and divide by weight in kg to get angstroms
         forceArr[i,j]=1E20*forceArr[i,j]*627.509*(4184/(0.529177*avNum))*(timestep^2)/(weight[i]/(avNum*1000))
         if ((diag>1) && (i==1)) print "i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3]" >> "diagnostics"
         if (diag>1) print i,weight[i],forceArr[i,1],forceArr[i,2],forceArr[i,3] >> "diagnostics"
         newarr[i,j]=oldarr[i,j]+damping*(oldarr[i,j]-olderarr[i,j])+forceArr[i,j]
         if ((i==fixedatom1) || (i==fixedatom2) || (i==fixedatom3) || (i==fixedatom4)) newarr[i,j]=oldarr[i,j]
#turn around atoms outside the box
         if (boxon==1) {
            if (newarr[i,j]>boxsize) if (oldarr[i,j]>olderarr[i,j]) newarr[i,j]=oldarr[i,j]+damping*(olderarr[i,j]-oldarr[i,j])+forceArr[i,j]
            if (newarr[i,j]<-1*boxsize) if (oldarr[i,j]<olderarr[i,j]) newarr[i,j]=oldarr[i,j]+damping*(olderarr[i,j]-oldarr[i,j])+forceArr[i,j]
            }
         }
      }
   }
########

if ((runpointnum % keepevery)==0) system("cat g09.log >> dyn")
print "%nproc=" processors
print "%mem=" memory
if (killcheck!=1) print "%chk=" checkpoint
if (nonstandard==0) {
   print "# " method " force scf=(tight,nosym) "
   if (meth2=="unrestricted") print "guess=mix" #for unrestricted calculations
   if (meth2=="read") print "guess=tcheck" #for reading orbitals from check, sometimes faster, sometimes not
   print "pop=none "
   if (length(meth3)>2) print meth3
   if (length(meth4)>2) print meth4
   }
if (nonstandard==1) {
   print "# "
   print "nonstd"
   system("cat nonstandard")
   }
print ""
print  title1,title2,title3,title4
print "runpoint ",runpointnum
print "runisomer ",isomernum
if (DRP==1) print "maxForce and forceMult and maxAtomMove",maxForce,forceMult,maxAtomMove
print ""
print charge,multiplicity
print numAtoms >> "traj"
print newPotentialE,title1,title2,title3,title4,"runpoint ",runpointnum,"runisomer ",isomernum >> "traj"
for (i=1;i<=numAtoms;i++) {
   printf("%s %.7f %.7f %.7f",atSym[i],newarr[i,1],newarr[i,2],newarr[i,3])
   printf("%s %.7f %.7f %.7f",atSym[i],newarr[i,1],newarr[i,2],newarr[i,3]) >> "traj"
   print "" >> "traj"
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
if ((nmrtype>0) && (nmrdo==1)) {
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
if ((nmrtype>1) && (nmrdo==1)) {
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
if ((nmrtype>2) && (nmrdo==1)) {
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

if ((nmrcc==1) && (nmrdo==1)) {
   print "CCSD(T) NMR calculation" > "ZMAT"
   for (i=1;i<=numAtoms;i++) {
      printf("%s %.7f %.7f %.7f",atSym[i],newarr[i,1],newarr[i,2],newarr[i,3]) >> "ZMAT"
      print "" >> "ZMAT"
      }
   print "" >> "ZMAT"
   print "*ACES2(CALC=CCSD[T],PROP=NMR,BASIS=dzp" >> "ZMAT"
   print "ABCDTYPE=AOBASIS,TREAT_PERT=SEQUENTIAL,CC_PROG=ECC" >> "ZMAT"
   print "COORD=CARTESIAN" >> "ZMAT"
   print "MEM_UNIT=GB,MEMORY=2)" >> "ZMAT"
   print "" >> "ZMAT"
   }  
}

