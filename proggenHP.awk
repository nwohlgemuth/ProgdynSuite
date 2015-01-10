BEGIN {
# aug 2013 summary of changes
#includes molecular rotation, ability to do multiple NMR calculations, ONIOM with link atoms, 
#nonstandard routes, handling of linear molecules using geometry linear, fixed but with atoms over 99 but
#bug varies with version of Gaussian, randomization based on PROCINFO (solved many problems), added initialDiss 3 for random
#phase of normal modes
# Aug 2010 changes classicalSpacing to 2 and upped possible excited states to 4000
# Jan 2009 - a number of little changes to improve reporting, precision, etc, specification of displacement on particular modes
# Jan 2009 cannonball trajectories.  adds desired energy to initial velocities based on file cannontraj, so one can shoot toward a ts
# updated Nov 2008 to incorporate running DRPs
# updated Nov 2008 to allow for start without an initial freq calc using classical = 2
# updated Aug 2008 added to atom list to handle a large number of atoms without changes needed
# updated June 2008 to incorportate new method for choosing displacements with initialdis 2
# updated Jan 17 2008 - bug fix for > 99 atoms, 300 excitations of low modes possible
# version August 2007 - incorporates classical trajectory calculation option
#also allows listing of number of imaginary frequencies
# version Sept 16, 2005 - incorportates searchdir but not yet rotation
# now reads random numbers from temp811, starting at a random place
# The input files are generated before this and are tempfreqs, tempredmass,
# tempfrc, tempmodes, and tempstangeos.
# It will count the number of atoms.  

# default parameters, including quassiclassical, no displacements, transition state, not a DRP
# do not change these - rather, change progdyn.conf to set the parameters
initialDis=0; timestep=1E-15; scaling=1.0; temp=298.15
classical=0; numimag=1; DRP=0; cannonball=0
charge=0; multiplicity=1; method="HF/3-21G"; memory=20000000
diag=1; checkpoint="g09.chk"; searchdir="positive"; boxon=0
boxsize=10; maxAtomMove=0.1; title1="you"; title2="need"
title3="a"; title4="progdyn.conf"; processors=1; highlevel=999
conver1=4.184E26 #dividing by this converts amu angs^2 /s^2 to kcal/mol
geometry="nonlinear";rotationmode=0

#initialization and constants
for (i=1;i<=10000;i++) {disMode[i]=-1}
i=1;j=1;k=1
c=29979245800; h=6.626075E-34; avNum=6.0221415E23
RgasK=0.00198588; RgasJ=8.31447
numAtoms=0; atomnumber=0; classicalSpacing=2
zpeGauss=0; zpeGaussK=0; zpePlusE=0; potentialE=0

# read progdyn.conf for configuration info
blankLineTester=10
while (blankLineTester>1) {
   getline < "progdyn.conf"
   if ($1=="method") method=$2
   if ($1=="charge") charge=$2
   if ($1=="multiplicity") multiplicity=$2
   if ($1=="memory") memory=$2
   if ($1=="processors") processors=$2
   if ($1=="checkpoint") checkpoint=$2
   if ($1=="diagnostics") diag=$2
   if ($1=="initialdis") initialDis=$2
   if ($1=="timestep") timestep=$2
   if ($1=="scaling") scaling=$2
   if ($1=="temperature") temp=$2
   if ($1=="searchdir") searchdir=$2
   if ($1=="classical") classical=$2
   if ($1=="numimag") numimag=$2
   if ($1=="geometry") geometry=$2
   if ($1=="highlevel") highlevel=$2
   if ($1=="boxon") boxon=$2
   if ($1=="boxsize") boxsize=$2
   if ($1=="DRP") DRP=$2; if (DRP==1) classical=2 #this lets one start a DRP from a point that is not a freq calc
   if ($1=="maxAtomMove") maxAtomMove=$2
   if ($1=="cannonball") cannonball=$2
   if ($1=="displacements") disMode[$2]=$3
   if ($1=="controlphase") controlPhase[$2]=$3
   if ($1=="rotationmode") rotationmode=$2
   if ($1=="title") {
      title1=$2
      title2=$3
      title3=$4
      title4=$5
      }
   blankLineTester=length($0)
   }

if (diag>=1) print "***************** starting proggen *****************" >> "diagnostics"
if (diag>=1) print "method,charge,multiplicity,memory" >> "diagnostics"
if (diag>=1) print method,charge,multiplicity,memory >> "diagnostics"
if (diag>=1) print "processors,checkpoint,title,initialdis,timestep,scaling,temperature" >> "diagnostics"
if (diag>=1) print processors,checkpoint,title1,title2,title3,title4,initialDis,timestep,scaling,temp >> "diagnostics"
if (diag>=1) print "classical,numimag,highlevel,boxon,boxsize,DRP,maxAtomMove,cannonball" >> "diagnostics"
if (diag>=1) print classical,numimag,highlevel,boxon,boxsize,DRP,maxAtomMove,cannonball >> "diagnostics"

# put geometries into array, also figure out number of atoms
# note that this picks out the last geometry in a file, assuming
# that if there is an optimization followed by a freq, nothing else follows
# kludgy - repeats last line twice - must be a better way
do {
   getline < "tempstangeos"
   if (oldline==$0) $0=""
   oldline=$0
   atom = $1
   if (atom>numAtoms) numAtoms=atom
   atNum[atom]=$2
   geoArr[atom,1]=$4; geoArr[atom,2]=$5; geoArr[atom,3]=$6
   geoArrOrig[atom,1]=$4; geoArrOrig[atom,2]=$5; geoArrOrig[atom,3]=$6
   velArr[atom,1]=0; velArr[atom,2]=0; velArr[atom,3]=0
   }
while (length($0) > 0)

#output the number of atoms, used in many routines
print numAtoms

# put in atomic symbols and atomic weights - assigns a default mass but then reads it from tempmasses when possible
for (i=1;i<=numAtoms;i++) {
   getline < "tempmasses"
   if (atNum[i]==1) {atSym[i]="H";atWeight[i]=1.00783}
   if (atNum[i]==2) {atSym[i]="He";atWeight[i]=4.0026}
   if (atNum[i]==3) {atSym[i]="Li";atWeight[i]=6.941}
   if (atNum[i]==4) {atSym[i]="Be";atWeight[i]=9.012}
   if (atNum[i]==5) {atSym[i]="B";atWeight[i]=10.811}
   if (atNum[i]==6) {atSym[i]="C";atWeight[i]=12.}
   if (atNum[i]==7) {atSym[i]="N";atWeight[i]=14.007}
   if (atNum[i]==8) {atSym[i]="O";atWeight[i]=15.9994}
   if (atNum[i]==9) {atSym[i]="F";atWeight[i]=18.9984}
   if (atNum[i]==10) {atSym[i]="Ne";atWeight[i]=20.1797}
   if (atNum[i]==11) {atSym[i]="Na";atWeight[i]=22.989}
   if (atNum[i]==12) {atSym[i]="Mg";atWeight[i]=24.305}
   if (atNum[i]==13) {atSym[i]="Al";atWeight[i]=26.98154}
   if (atNum[i]==14) {atSym[i]="Si";atWeight[i]=28.0855}
   if (atNum[i]==15) {atSym[i]="P";atWeight[i]=30.9738}
   if (atNum[i]==16) {atSym[i]="S";atWeight[i]=32.066}
   if (atNum[i]==17) {atSym[i]="Cl";atWeight[i]=35.4527}
   if (atNum[i]==18) {atSym[i]="Ar";atWeight[i]=39.948}
   if (atNum[i]==19) {atSym[i]="K";atWeight[i]=39.0983}
   if (atNum[i]==20) {atSym[i]="Ca";atWeight[i]=40.078}
   if (atNum[i]==21) {atSym[i]="Sc";atWeight[i]=44.96}
   if (atNum[i]==22) {atSym[i]="Ti";atWeight[i]=47.867}
   if (atNum[i]==23) {atSym[i]="V";atWeight[i]=50.94}
   if (atNum[i]==24) {atSym[i]="Cr";atWeight[i]=51.9961}
   if (atNum[i]==25) {atSym[i]="Mn";atWeight[i]=54.938}
   if (atNum[i]==26) {atSym[i]="Fe";atWeight[i]=55.845}
   if (atNum[i]==27) {atSym[i]="Co";atWeight[i]=58.933}
   if (atNum[i]==28) {atSym[i]="Ni";atWeight[i]=58.693}
   if (atNum[i]==29) {atSym[i]="Cu";atWeight[i]=63.546}
   if (atNum[i]==30) {atSym[i]="Zn";atWeight[i]=65.38}
   if (atNum[i]==31) {atSym[i]="Ga";atWeight[i]=69.723}
   if (atNum[i]==32) {atSym[i]="Ge";atWeight[i]=72.64}
   if (atNum[i]==33) {atSym[i]="As";atWeight[i]=74.9216}
   if (atNum[i]==34) {atSym[i]="Se";atWeight[i]=78.96}
   if (atNum[i]==35) {atSym[i]="Br";atWeight[i]=79.904}
   if (atNum[i]==46) {atSym[i]="Pd";atWeight[i]=106.42}
   if (atNum[i]==53) {atSym[i]="I";atWeight[i]=126.90447}
# gets actual weight from freqinHP when possible so a prior calc with readisotopes gets you isotopic substitution
   if ((i<100) && ($9>0)) atWeight[i]=$9
#  if ((i>99) && ($8>0)) atWeight[i]=$8

   if ((diag>1) && (i==1)) print "atNum[i],atSym[i],atWeight[i],geoArr[i,1],geoArr[i,2],geoArr[i,3]" >> "diagnostics"
   if (diag>1) print atNum[i],atSym[i],atWeight[i],geoArr[i,1],geoArr[i,2],geoArr[i,3] >> "diagnostics"
   }

# read in frequencies, scale them, read in Reduced masses, read in force 
#constants, replace negative frequencies by 2 wavenumbers
numFreq=3*numAtoms-6
if (geometry=="linear") numFreq=3*numAtoms-5
for (i=1;i<=numFreq;i++) {
   $0=""
   getline < "tempfreqs"
   freq[i]=$0*scaling
   if (freq[i]<0) freq[i]=2
   }
for (i=1;i<=numFreq;i++) {
   $0=""
   getline < "tempredmass"
   redMass[i]=$0
   if (redMass[i]=="") redMass[i]=1.
   }
for (i=1;i<=numFreq;i++) {
   $0=""
   getline < "tempfrc"
   frc[i]=$0
   if (frc[i]=="") frc[i]=0.0001
   if (frc[i]==0) frc[i]=0.0001
   if ((diag>1) && (i==1)) print "freq[i],redMass[i],frc[i]" >> "diagnostics"
   if (diag>1) print freq[i],redMass[i],frc[i] >> "diagnostics"
   }

# read in the modes - note that trajectories always need a freq calc with freq=hpmodes unless classical=2
if (classical!=2) {
   for (i=1;i<=numFreq;i+=5) {
      for (j=1;j<=(3*numAtoms);j++) {
         getline < "tempmodes"
         mode[i,$2,$1]=$4; mode[i+1,$2,$1]=$5; mode[i+2,$2,$1]=$6; mode[i+3,$2,$1]=$7; mode[i+4,$2,$1]=$8
         }
      }
   }
if (diag>2) {for (i=1;i<=numFreq;i++) {print mode[i,1,1],mode[i,1,2],mode[i,1,3] >> "modesread"}}

# if doing a cannonball trajectory, read in the vector
if (cannonball>0) {
   for (i=1;i<=numAtoms;i++) {
      getline < "cannontraj"
      cannonArr[i,1]=$1; cannonArr[i,2]=$2; cannonArr[i,3]=$3
      }
   }

# collect a series of random numbers from file temp811, generated from an outside random number generator called by prodynstarterHP
# read from temp811, starting at a random place
srand(PROCINFO["pid"]); tester=rand()*1000
for (i=1;i<=tester;i++) getline < "temp811"
for (i=1;i<=numFreq;i++) {
   getline < "temp811"; randArr[i]=$1
   getline < "temp811"; randArrB[i]=$1
   getline < "temp811"; randArrC[i]=$1
}
for (i=1;i<=6;i++) {
   getline < "temp811"; randArrR[i]=$1
}

# for a QM distribution for a harmonic oscillator in its ground state, we want to generate a set of random numbers
#between -1 and 1 weighted such that numbers toward the center are properly more common
i=1
while (i<=numFreq) {
   if ((initialDis==2) || (disMode[i]==2)) {
      getline < "temp811"
      tempNum=2*($1-.5)
      prob=exp(-(tempNum^2))
      getline < "temp811"
      if ($1<prob) {
         randArrD[i]=tempNum
         i++
         }
      }
   if ((initialDis!=2) && (disMode[i]!=2)) i++
   }

# to start without normal modes or frequencies we need to just pick a random direction for the motion of each atom, requiring 3N random numbers
for (i=1;i<=numAtoms;i++) {
   for (j=1;j<=3;j++) {
   getline < "temp811"
   if ($1>0.5) randArrE[i,j]=1
   if ($1<.5) randArrE[i,j]=-1
      }
   }

# determine energy in each normal mode
for (i=1;i<=numFreq;i++) {
   zpeJ[i]=0.5*h*c*freq[i]       #units J per molecule
#if classical, treat as modes spaced by classicalSpacing wavenumbers
   if (classical==1) zpeJ[i]=0.5*h*c*classicalSpacing  # the zpe is not used when classical but the spacing is used to calculate the E in mode
   zpeK[i]=zpeJ[i]*avNum/4184    #units kcal/mol
   if (temp<10) vibN[i]=0        # avoids working with very small temperatures - if the temp is too low, it just acts like 0 K
   if (temp>=10) {
      zpeRat[i]=exp((-2*zpeK[i])/(RgasK*temp))
      if (zpeRat[i]==1) zpeRat[i]=.99999999999
      Q[i]=1/(1-zpeRat[i])
      newRand=randArr[i]
      vibN[i]=0
      tester=1/Q[i]
#     get up to 4000 excitations of low modes
      for (j=1;j<=(4000*zpeRat[i]+2);j++) {
         if (newRand>tester) vibN[i]++
         tester=tester+((zpeRat[i]^j)/Q[i])
         }
      }
   }

# figure out mode energies and maximum classical shift and then actual shift
# also calculated total energy desired for molecule
desiredModeEnK=0
for (i=1;i<=numFreq;i++) {
   modeEn[i]=(zpeJ[i]*1E18)*(2*vibN[i]+1) # units here are mDyne Ansgroms for compatability with Gaussian force constants
   if (classical==1) modeEn[i]=(zpeJ[i]*1E18)*2*vibN[i]    #no zpe when classical
   modeEnK[i]=zpeK[i]*(2*vibN[i]+1)
   if (classical==1) modeEnK[i]=zpeK[i]*2*vibN[i]          #no zpe when classical
   desiredModeEnK=desiredModeEnK + modeEnK[i]
# no 1/2 hv for imaginary frequencies
# treating modes with frequencies <10 as translations, ignoring their zero point energies
   if (freq[i]<10) modeEn[i]=(zpeJ[i]*1E18)*(2*vibN[i])
   maxShift[i]=(2*modeEn[i]/frc[i])^0.5
# new 2012 initialDis 3 means random phase of normal mode
   if (initialDis==3) shift[i]=maxShift[i]*sin(randArrC[i]*3.141592*2)
   if (initialDis==2) shift[i]=maxShift[i]*randArrD[i]
   if (initialDis==1) shift[i]=maxShift[i]*(2*(randArrC[i]-0.5))
   if (initialDis==0) shift[i]=0
# lines below allow for setting of displacement mode for individual modes
# It used to be necessary to use disMode 10 to turn off displacements for a mode, but hopefully that bug is killed and you can use disMode 0
   if (disMode[i]==3) shift[i]=maxShift[i]*sin(randArrC[i]*3.141592*2)
   if (disMode[i]==2) shift[i]=maxShift[i]*randArrD[i]
   if (disMode[i]==1) shift[i]=maxShift[i]*(2*(randArrC[i]-0.5))
   if (disMode[i]==10) shift[i]=0 #kept for backward compatability
   if (disMode[i]==0) shift[i]=0
# no displacements along imaginary frequencies and very low ones - it is better to treat these
# as translations - employing a shift can give you initial weird geometries
   if (freq[i]<10) shift[i]=0
   if (numimag==1) shift[1]=0
   if (numimag==2) shift[2]=0
   }
for (i=1;i<=numFreq;i++) {
   if ((diag>1) && (i==1)) print "zpeJ[i],zpeK[i],zpeRat[i],Q[i],vibN[i],modeEn[i],maxShift[i],shift[i]" >> "diagnostics"
   if (diag>1) print zpeJ[i],zpeK[i],zpeRat[i],Q[i],vibN[i],modeEn[i],maxShift[i],shift[i] >> "diagnostics"
   }

# multiply each of the modes by its shift and add them up
# Do not do this if classical=2
if (classical!=2) {
   for (i=1;i<=numFreq;i++) {
      for (j=1;j<=numAtoms;j++) {
         for (k=1;k<=3;k++) {
            shiftMode[i,j,k]=mode[i,j,k]*shift[i]
            geoArr[j,k]=geoArr[j,k]+shiftMode[i,j,k]
            }
         }
      }
   }

#now start toward velocities
for (i=1;i<=numFreq;i++) {
   kinEn[i]=100000*(modeEn[i]-0.5*frc[i]*shift[i]^2)  # the 100000 converts to g angstrom^2 s^2
   vel[i]=(2*kinEn[i]/(redMass[i]/avNum))^0.5        # in angstrom / s
#use searchdir in progdyn.conf to control the direction for trajectories started from a saddle point
   if (numimag>1) numimag=1  #only the first freq can be sent in the searchdir direction, the rest go in a random direction
   if (i>numimag) {
      if (randArrB[i]<0.5) vel[i]=-vel[i]
      }
   if (i==numimag) {
      if (searchdir=="negative") vel[i]=-vel[i]
      }
   if ((diag>1) && (i==1)) print "vel[i]" >> "diagnostics"
   if (diag>1) print vel[i] >> "diagnostics"
   }

# if controlphase is being used, set the velocity on particular modes as positive or negative as requested
for (i=1;i<=numFreq;i++) {
   if ((controlPhase[i]=="positive") && (vel[i]<0)) vel[i]=-vel[i]
   if ((controlPhase[i]=="negative") && (vel[i]>0)) vel[i]=-vel[i]
   }

# multiply each of the modes by its velocity and add them up
# Do not do this if classical=2
if (classical!=2) {
   for (i=1;i<=numFreq;i++) {
      for (j=1;j<=numAtoms;j++) {
         for (k=1;k<=3;k++) {
            velMode[i,j,k]=mode[i,j,k]*vel[i]*timestep
            velArr[j,k]=velArr[j,k]+velMode[i,j,k]
            }
         }
      }
   }

# to start without normal modes or frequencies we figure out the energy per atom based on 1/2RT in degree of freedom
if (classical==2) {
# to avoid a bug with a box on, starts without modes should use the input geometry, not the standard
   do {
      getline < "tempinputgeos"
      if (oldline==$0) $0=""
      oldline=$0
      atom = $1
      geoArr[atom,1]=$4; geoArr[atom,2]=$5; geoArr[atom,3]=$6
      geoArrOrig[atom,1]=$4; geoArrOrig[atom,2]=$5; geoArrOrig[atom,3]=$6
      }
   while (length($0) > 0)
   degFreedomEnK=temp*RgasK
   degFreedomEnJ=degFreedomEnK/(avNum/4184)
   cartEn=degFreedomEnJ*1E18
   kinEnCart=100000*cartEn
#print degFreedomEnK, degFreedomEnJ, cartEn, kinEnCart
   for (i=1;i<=numAtoms;i++) {
      for (j=1;j<=3;j++) {
         velArr[i,j]=randArrE[i,j]*timestep*(2*kinEnCart/(atWeight[i]/avNum))^0.5
         if (DRP==1) velArr[i,j]=0
         }
      }
   }

# calculate the KE in the modes at this point
KEinitmodes=0
for (j=1;j<=numAtoms;j++) {
   KEinitmodes=KEinitmodes + 0.5*atWeight[j]*(velArr[j,1]^2 + velArr[j,2]^2 + velArr[j,3]^2)/((timestep^2)*conver1)
   }

# add molecular rotation if requested
if (rotationmode>0) {
#establish three rotation vectors
   for (j=1;j<=numAtoms;j++) {
      rotateX[j,1]=0
      rotateX[j,2]=-geoArrOrig[j,3]
      rotateX[j,3]=geoArrOrig[j,2]
      rotateY[j,1]=-geoArrOrig[j,3]
      rotateY[j,2]=0
      rotateY[j,3]=geoArrOrig[j,1]
      rotateZ[j,1]=-geoArrOrig[j,2]
      rotateZ[j,2]=geoArrOrig[j,1]
      rotateZ[j,3]=0
      }
#figure out how much energy is in the raw vectors
   eRotX=0;eRotY=0;eRotZ=0
   for (j=1;j<=numAtoms;j++) {
      for (k=1;k<=3;k++) {
         eRotX=eRotX + 0.5*atWeight[j]*(rotateX[j,k]^2)/((timestep^2)*conver1)
         eRotY=eRotY + 0.5*atWeight[j]*(rotateY[j,k]^2)/((timestep^2)*conver1)
         eRotZ=eRotZ + 0.5*atWeight[j]*(rotateZ[j,k]^2)/((timestep^2)*conver1)
         }
      }
#  print "rotation energies if raw vector used",eRotX,eRotY,eRotZ
#now deciie how much energy we want in each rotation
   keRx=-0.5*0.001987*temp*log(1-randArrR[1])
   keRy=-0.5*0.001987*temp*log(1-randArrR[2])
   keRz=-0.5*0.001987*temp*log(1-randArrR[3])
   if (eRotX<1) keRx=0;if (eRotY<1) keRy=0;if (eRotZ<1) keRz=0
   rotEdesired=keRx+keRy+keRz
   signX=1;signY=1;signZ=1
   if (randArrR[4]<.5) signX=-1
   if (randArrR[5]<.5) signY=-1
   if (randArrR[6]<.5) signZ=-1

#  print "desired energies",keRx,keRy,keRz,"and random numbers",randArrR[1],randArrR[2],randArrR[3]
#protect against zero rotations 
   if (eRotX<1) eRotX=1;if (eRotY<1) eRotY=1;if (eRotZ<1) eRotZ=1
#now scale the rotational vectors
   scaleX=(keRx/eRotX)^.5
   scaleY=(keRy/eRotY)^.5
   scaleZ=(keRz/eRotZ)^.5
#  print "scaling factors" scaleX,scaleY,scaleZ
   for (j=1;j<=numAtoms;j++) {
      for (k=1;k<=3;k++) {
         rotateX[j,k]=rotateX[j,k]*scaleX*signX
         rotateY[j,k]=rotateY[j,k]*scaleY*signY
         rotateZ[j,k]=rotateZ[j,k]*scaleZ*signZ
         } 
      }
   for (j=1;j<=numAtoms;j++) {
#     print rotateX[j,1],"   ",rotateX[j,2],"   ",rotateX[j,3]
      }
#  print ""
   for (j=1;j<=numAtoms;j++) {
#     print rotateY[j,1],"   ",rotateY[j,2],"   ",rotateY[j,3]
      }
#  print ""
   for (j=1;j<=numAtoms;j++) {
#     print rotateZ[j,1],"   ",rotateZ[j,2],"   ",rotateZ[j,3]
      }
# now add the rotational vectors  to velArr
   for (j=1;j<=numAtoms;j++) {
      for (k=1;k<=3;k++) {
         velArr[j,k]=velArr[j,k]+rotateX[j,k]+rotateY[j,k]+rotateZ[j,k]
         }
      }
   }

# if doing a cannonball, adjust multiplier until extra energy is correct
if (cannonball>0) {
   multiplier=1; tester=0; tolerance=.1
   while (tester==0) {
      KEinittotal=0
      for (j=1;j<=numAtoms;j++) {
         cannonvelArr[j,1]=velArr[j,1]+multiplier*cannonArr[j,1]; cannonvelArr[j,2]=velArr[j,2]+multiplier*cannonArr[j,2]; cannonvelArr[j,3]=velArr[j,3]+multiplier*cannonArr[j,3]
         KEinittotal=KEinittotal + 0.5*atWeight[j]*(cannonvelArr[j,1]^2 + cannonvelArr[j,2]^2 + cannonvelArr[j,3]^2)/((timestep^2)*conver1)
         }
      if (KEinittotal>(KEinitmodes+cannonball+tolerance)) multiplier=multiplier*0.98901364
      if (KEinittotal<(KEinitmodes+cannonball-tolerance)) multiplier=multiplier*1.01
      if ((KEinittotal<(KEinitmodes+cannonball+tolerance)) && (KEinittotal>(KEinitmodes+cannonball-tolerance))) tester=1
      }
   for (j=1;j<=numAtoms;j++) {
      velArr[j,1]=velArr[j,1]+multiplier*cannonArr[j,1]; velArr[j,2]=velArr[j,2]+multiplier*cannonArr[j,2]; velArr[j,3]=velArr[j,3]+multiplier*cannonArr[j,3]
      }
   }

#output the new geometry.
# ****** this section changed for special experiment for cyclopentadiene.  do not use this for other cases
# atWeight[4]=140.0001
# ****** line below added for special experiment switching mass from 12 to 140, keeping momenta the same
#velArr[4,1]=velArr[4,1]/11.66667; velArr[4,2]=velArr[4,2]/11.66667; velArr[4,3]=velArr[4,3]/11.66667
for (j=1;j<=numAtoms;j++) {
   printf("%2s % .7f % .7f % .7f %9.5f \n",atSym[j],geoArr[j,1],geoArr[j,2],geoArr[j,3],atWeight[j])
   }

#output the velocities and calculate the total kinetic energy overall
KEinittotal=0
for (j=1;j<=numAtoms;j++) {
   KEinittotal=KEinittotal + 0.5*atWeight[j]*(velArr[j,1]^2 + velArr[j,2]^2 + velArr[j,3]^2)/((timestep^2)*conver1)
   printf("% .8f % .8f % .8f \n",velArr[j,1],velArr[j,2],velArr[j,3])
   }

#anything else I add to the file will not affect the trajectories but will keep a record and be good for analysis
for (i=1;i<=numFreq;i++) {
   if (initialDis==0) printf("%.6f   % .6f    %4i    % 1.4e       % .6f %1i\n", randArr[i], randArrB[i], vibN[i], vel[i], shift[i], disMode[i])
   if (initialDis==1) printf("%.6f   % .6f    %4i    % 1.4e       % .6f %1i\n", randArr[i], randArrC[i], vibN[i], vel[i], shift[i], disMode[i])
   if (initialDis==2) printf("%.6f   % .6f    %4i    % 1.4e       % .6f %1i\n", randArr[i], randArrD[i], vibN[i], vel[i], shift[i], disMode[i])
   if (initialDis==3) printf("%.6f   % .6f    %4i    % 1.4e       % .6f %1i % .6f\n", randArr[i], randArrC[i], vibN[i], vel[i], shift[i], disMode[i], sin(randArrC[i]*3.141592*2))
   }
print "temp ",temp
print "initialDis",initialDis
print "classical",classical
print "timestep",timestep
print "numimag",numimag
OFMT = "%.3f"
print "Total mode energy desired=",desiredModeEnK
print "KE initial from modes=",KEinitmodes,"   KE initial total=",KEinittotal,"   Rotational Energy desired=",rotEdesired
if (cannonball>0) print "cannonball",cannonball,"  cannon Energy=",KEinittotal-KEinitmodes
if (boxon>0) print "boxsize",boxsize
if (DRP>0) print "DRP",DRP,"   maxAtomMove",maxAtomMove
if (DRP>0) print maxAtomMove > "maxMove"
}  # End of BEGIN

/Zero-point correction/ {zpeGauss=$3}
/zero-point Energies/ {zpePlusE=$7}
END {
zpeGaussK=zpeGauss*627.509
potentialE=zpePlusE - zpeGauss
OFMT = "%.6f"
print "Gaussian zpe=",zpeGauss,"or",zpeGaussK,"kcal/mol  E + zpe=",zpePlusE,"  potential E=",potentialE
print "" #will use blank line to mark end of geoPlusVel file
}

