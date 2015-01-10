# this program must be modified for each system.  It interprets the gaussian output and writes
# data to dynfollowfile
BEGIN {
firsttitle=1
getline < "isomernumber"
isomer=$1
}
# The word here should be the first word of the title.
/ alame3/ {
   if (firsttitle==1) {
      printf("%s %s %s %s %s %s %s  ",$1,$2,$3,$4,$6,$7,$8)
      runpoint=$6
      }
   firsttitle++
   }
/Standard orientation/,/Rotational constants/ {
   if (($1>.5) && ($1<30)) {
      A[$1]=$4;B[$1]=$5;C[$1]=$6
      }
   }
/before annihilation/ {
   printf("%s %.5f ",$1,$6)
   }

END {
#vary distances checked and response by system
   C1N6=Distance(1,6)
   C3C8=Distance(3,8)
   C1C8=Distance(1,8)
   printf("%s %.3f  %s %.3f  %s %.3f ","C1N6",C1N6,"C3C8",C3C8,"C1C8",C1C8)
   if (runpoint>500) {
      print "Too many points.  XXXXN"
#     system("date > nogo")
      }
      if ((C1N6>2.5) && (C3C8<1.8)) {
         getline < "skipstart"
#        if ($1=="reverse") system("date > nogo")
         print "2,3-Rearrangement XXXX23N"
         }
      if ((C1N6>2.5) && (C1C8<1.8)) {
         getline < "skipstart"
#        if ($1=="reverse") system("date > nogo")
         print "1,2-Rearrangement XXXX12N"
         }
      if ((C1N6>3.2) && (C3C8>3.4) && (C1C8>3.2)) {
         getline < "skipstart"
#        if ($1=="reverse") system("date > nogo")
         print "Dissociation XXXXDN"
         }
      if ((C1N6<1.7) && (C3C8>2.8)) {
         getline < "skipstart"
#        if ($1=="reverse") system("date > nogo")
         print "reformed starting material XXXX23N"
         }
   system("date '+%b:%d:%Y %T'")
   system("tail -1 Echeck | grep XXXX")
   }

#the rest of the file is just functions to help in checking the geometry
function Distance(Atom1,Atom2) {
  return sqrt((A[Atom1]-A[Atom2])^2 + (B[Atom1]-B[Atom2])^2 + (C[Atom1]-C[Atom2])^2)
}

function Angle(Atom1,Atom2,Atom3) {
   value=((-Distance(Atom1,Atom3)^2+Distance(Atom1,Atom2)^2+Distance(Atom2,Atom3)^2)/(2*Distance(Atom1,Atom2)*Distance(Atom2,Atom3)))
   return acos(value)
}

function asin(x) { return (180/3.141592)*atan2(x, sqrt(1-x*x)) }

function acos(x) { return (180/3.141592)*atan2(sqrt(1-x*x), x) }

function atan(x) { return (180/3.141592)*atan2(x,1) }

function Dihedral(Atom1,Atom2,Atom3,Atom4) {
   B1x=A[Atom2]-A[Atom1]
   B1y=B[Atom2]-B[Atom1]
   B1z=C[Atom2]-C[Atom1]
   B2x=A[Atom3]-A[Atom2]
   B2y=B[Atom3]-B[Atom2]
   B2z=C[Atom3]-C[Atom2]
   B3x=A[Atom4]-A[Atom3]
   B3y=B[Atom4]-B[Atom3]
   B3z=C[Atom4]-C[Atom3]
   modB2=sqrt((B2x^2)+(B2y^2)+(B2z^2))
# yAx is x-coord. etc of modulus of B2 times B1
   yAx=modB2*(B1x)
   yAy=modB2*(B1y)
   yAz=modB2*(B1z)
# CP2 is the crossproduct of B2 and B3
   CP2x=(B2y*B3z)-(B2z*B3y)
   CP2y=(B2z*B3x)-(B2x*B3z)
   CP2z=(B2x*B3y)-(B2y*B3x)
   termY=((yAx*CP2x)+(yAy*CP2y)+(yAz*CP2z))
# CP is the crossproduct of B1 and B2
   CPx=(B1y*B2z)-(B1z*B2y)
   CPy=(B1z*B2x)-(B1x*B2z)
   CPz=(B1x*B2y)-(B1y*B2x)
   termX=((CPx*CP2x)+(CPy*CP2y)+(CPz*CP2z))
  dihed4=(180/3.141592)*atan2(termY,termX)
  return dihed4
}

function killdyn(isomer) {
   system("rm -f dyn")
}


