#!/bin/bash
#progdynstarterHP, made to use high-precision modes from Gaussian output with freq=hpmodes
#updated to create a random number file temp811 that is used by proggenHP
#version Sep 2005, made for workstations
#version Aug 2007 to allow periodic copying of g09.log to dyn putting it under control of progdynb
#version Feb 2008 moves variables like the scratch directory and location of randgen to the beginning
#version Mar 2008 added proganal reporting to points 1 and 2
#version Jan 2009 fixed bug generator of having proganal run twice in checking for complete runs
#version May 2009 Echeck catches bad energies after only one point, other lines written simpler, triple while loop, revised comments
#version Aug 2010 isomernumber adds words to ease parsing, increased elements up to bromine, runpointnumber checked for more appropriate restarts
#version Aug 2011 runpointnumber starts better, restart better if died during first few points, awk bug fix
#version Aug 2012 freqinHP reads with only 3 freqs, goingwell and other temp files moved to $SCRATC_DIR
#version Aug 2013 adds ability to automatically run a CFOUR program if the file ZMAT exists
#
#version Oct 2014 updated to run correctly on BYU's FSL 
#LIMITATIONS - standard version only handles elements up to bromine, must change program to do higher atomic numbers
#   only handles up to 4000th excited state for modes - this could start to affect the initialization of classical modes or transition vectors at
#    extremely high temperatures
#   The routine that checks whether the actual energy approximately equals the desired energy checks for lines containing "SCF Done" or "EUMP2 =" or " Energy="
#   This should handle ordinary calculations HF, DFT, ONIOM, and MP2 calculatons but the routine in prog2ndpoint would have to be changed for other calcs. 
#
#                                        OUTLINE
# A. initilize to perform Gaussian jobs and know where we are
#    start loop
# B. if no file named "skipstart" then generate a new isomer.  Instructions: Get rid of skipstart to start new isomer.
#    the B loop generates geoPlusVel, adds it to geoRecord, generates and runs first and second points, and sets up for continuous loop
# C. loop over propagation steps
# 
#  AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#TEMPORARY_DIR, RAND_DIR, scratchdir, G09_ROOT, LOG_FILE all may need varied from system to system and assigned here or by program calling this one

#Variable Setup
LOG_FILE="docslog"
RUN_NUMBER=0

# File System Setp
export LC_ALL=C
export JOB_NAME=freqinHP
#export TEMPORARY_DIR=
export GAUSS_SCRDIR="$TEMPORARY_DIR/temporary_files"
export PROG_SCRDIR="$TEMPORARY_DIR/prog_files"
export PROG_HOME="$TEMPORARY_DIR"
export G09_ROOT="/blueapps/chem"
export RAND_DIR="$TEMPORARY_DIR"
#. $G09_ROOT/gaussian/g09/bsd/g09.profile

cd $TEMPORARY_DIR

rm -f nogo    # assume that if someone is starting a job, they want it to go.
rm -f diagnostics goingwell tempdone # diagnostics contains extra info from previous runs, other two files are from older versions of progdyn

#### Triple 'while' loop - will have to break multiple times to get out, but advantage is ability to control starting over
while (true)
do

# As long as there is a file "goingwell" the program will not exit entirely by itself
rm -f $PROG_SCRDIR/goingwell
while (true)
do
#  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
   if (test -f skipstart) then
      echo "skipping start and continuing from previous runs"
   else
#  B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1B1 generate geoPlusVel and first input file
      if [ "`cat runpointnumber`" = "1" ]; then
         echo "X did not complete first point so new isomer started" >> dynfollowfile
      fi
      if [ "`cat runpointnumber`" = "2" ]; then
         echo "X did not complete second point so new isomer started" >> dynfollowfile
      fi
      if [ "`cat runpointnumber`" = "3" ]; then
         echo "X did not complete third point so new isomer started" >> dynfollowfile
      fi
      
      if (test -f bypassproggen) then
         echo "taking starting conditions from pre-generated geoPlusVel"
      else
         $RAND_DIR/randgen > temp811
# the next 8 lines would have to be changed to use low-precision modes
         awk '/        1         2/,/Harmonic frequencies/ {print}' $JOB_NAME > temp401
         awk '/Frequencies --/ {print $3;print $4;print $5;print $6;print $7}' temp401 > tempfreqs
         awk '/Reduced masses/ {print $4;print $5;print $6;print $7;print $8}' temp401 > tempredmass
         awk '/Force constants/ {print $4;print $5;print $6;print $7;print $8}' temp401 > tempfrc
         awk '/0/ && ((length($1) < 2) && ($1 < 4)) {print}' temp401 > tempmodes
         awk '/has atomic number/ {print}' $JOB_NAME > tempmasses
         awk '/Standard orientation:/,/tional const/ {if ($3==0) print}' $JOB_NAME > tempstangeos
         awk '/Input orientation:/,/Stoichiometry/ {if ($3==0) print}' $JOB_NAME > tempinputgeos
         awk -f $TEMPORARY_DIR/proggenHP $JOB_NAME > geoPlusVel
      fi
      if (test -f isomernumber) then
         cp isomernumber temp533
         awk 'BEGIN {getline;i=$1+1;print i,"----trajectory isomer number----"}' temp533 > isomernumber
         rm temp533
      else
         echo "1 ----trajectory isomer number----" > isomernumber
      fi
      echo 1 > runpointnumber
      if (test -f g09.com) then
        rm g09.com
      fi
      awk -f $TEMPORARY_DIR/prog1stpoint isomernumber > g09.com
#  B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2B2  if first part successfule then clean up and run the first input file, otherwise die
      if (test -s g09.com) then
         rm tempfreqs tempredmass tempfrc tempmodes tempstangeos tempmasses temp401 temp811 tempinputgeos
         cat isomernumber >> geoRecord
         cat geoPlusVel >> geoRecord
         rm -f $PROG_SCRDIR/goingwell
         cp $TEMPORARY_DIR/g09.com $PROG_SCRDIR/g09.com
	 cp $PROG_SCRDIR/g09.com $COM_LOG_FILES/g09"$RUN_NUMBER".com
	 $G09_ROOT/bin/rung09 $PROG_SCRDIR/g09.com > $PROG_SCRDIR/g09.log
	 cp $PROG_SCRDIR/g09.log $COM_LOG_FILES/g09"$RUN_NUMBER".log
	 ((RUN_NUMBER++))
         grep 'Normal termination' $PROG_SCRDIR/g09.log > $PROG_SCRDIR/goingwell
         if (test -s $PROG_SCRDIR/goingwell) then
            cat $PROG_SCRDIR/g09.log >> dyn
            cp $PROG_SCRDIR/g09.log olderdynrun
         else
            cp $PROG_SCRDIR/g09.log $TEMPORARY_DIR/g09.log
            break
         fi
      else
         break
      fi
#  B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3B3 if B2 worked then you are here.  create 2nd point, run it, and set up for propagation loop
      rm g09.com
      echo 2 > runpointnumber
      awk -f $TEMPORARY_DIR/prog2ndpoint $PROG_SCRDIR/g09.log > g09.com
# before we decide to run this, check the energy
      awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> dynfollowfile
      rm -f $PROG_SCRDIR/tempdone
      tail -1 dynfollowfile | awk '/XXXX/ {print}' > $PROG_SCRDIR/tempdone
      if (test  -s $PROG_SCRDIR/tempdone) then
         rm -f dyn
         rm -f traj
         echo 0 > runpointnumber
         break
      fi
      if (test -s g09.com) then
         rm -f $PROG_SCRDIR/goingwell
         cp $TEMPORARY_DIR/g09.com $PROG_SCRDIR/g09.com
	 cp $PROG_SCRDIR/g09.com $COM_LOG_FILES/g09"$RUN_NUMBER".com
	 $G09_ROOT/bin/rung09 $PROG_SCRDIR/g09.com > $PROG_SCRDIR/g09.log
	 cp $PROG_SCRDIR/g09.log $COM_LOG_FILES/g09"$RUN_NUMBER".log
	 ((RUN_NUMBER++))
         grep 'Normal termination' $PROG_SCRDIR/g09.log > $PROG_SCRDIR/goingwell
         if (test -s $PROG_SCRDIR/goingwell) then
            cp $PROG_SCRDIR/g09.log olddynrun
            cat $PROG_SCRDIR/g09.log >> dyn
            awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> dynfollowfile
            awk '/Input orientation/,/Distance matrix/ {print}' olddynrun | awk '/   0   / {print}' > old
            awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
            echo 3 > runpointnumber
            awk -f $TEMPORARY_DIR/progdynb olddynrun > g09.com
            rm -f old older 
         else
            cp $PROG_SCRDIR/g09.log $TEMPORARY_DIR/g09.log
            break
         fi
      else
         break
      fi
# we've just completed a start, so lets skipstart until instructed otherwise
      echo "forward" > skipstart
   fi
# Reverse trajectories starter routine
   if [ `cat skipstart` = "reverserestart" ]; then
      
      rm g09.com
      echo 1 > runpointnumber
      awk -f $TEMPORARY_DIR/prog1stpoint isomernumber > g09.com
      if (test -s g09.com) then
         rm -f $PROG_SCRDIR/goingwell
         cp $TEMPORARY_DIR/g09.com $PROG_SCRDIR/g09.com
	 cp $PROG_SCRDIR/g09.com $COM_LOG_FILES/g09"$RUN_NUMBER".com
	 $G09_ROOT/bin/rung09 $PROG_SCRDIR/g09.com > $PROG_SCRDIR/g09.log
	 cp $PROG_SCRDIR/g09.log $COM_LOG_FILES/g09"$RUN_NUMBER".log
	 ((RUN_NUMBER++))
         grep 'Normal termination' $PROG_SCRDIR/g09.log > $PROG_SCRDIR/goingwell
         if (test -s $PROG_SCRDIR/goingwell) then
            cp $PROG_SCRDIR/g09.log olderdynrun
         else
            cp $PROG_SCRDIR/g09.log $TEMPORARY_DIR/g09.log
            break
         fi
      else
         break
      fi
      rm g09.com
      echo 2 > runpointnumber
      awk -f $TEMPORARY_DIR/prog2ndpoint $PROG_SCRDIR/g09.log > g09.com
      awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> dynfollowfile
      rm -f $PROG_SCRDIR/tempdone
      if (test -s g09.com) then
         rm -f $PROG_SCRDIR/goingwell
         cp $TEMPORARY_DIR/g09.com $PROG_SCRDIR/g09.com
         cp $PROG_SCRDIR/g09.com $COM_LOG_FILES/g09"$RUN_NUMBER".com
	 $G09_ROOT/bin/rung09 $PROG_SCRDIR/g09.com > $PROG_SCRDIR/g09.log
	 cp $PROG_SCRDIR/g09.log $COM_LOG_FILES/g09"$RUN_NUMBER".log
	 ((RUN_NUMBER++))
         grep 'Normal termination' $PROG_SCRDIR/g09.log > $PROG_SCRDIR/goingwell
         if (test -s $PROG_SCRDIR/goingwell) then
            cp $PROG_SCRDIR/g09.log olddynrun
            cat $PROG_SCRDIR/g09.log >> dyn
            awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> dynfollowfile
            awk '/Input orientation/,/Distance matrix/ {print}' olddynrun | awk '/   0   / {print}' > old
            awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
            echo 3 > runpointnumber
            awk -f $TEMPORARY_DIR/progdynb olddynrun > g09.com
            rm -f old older
         else
            cp $PROG_SCRDIR/g09.log $TEMPORARY_DIR/g09.log
            break
         fi
      else
         break
      fi
# we've just completed a reversestart, so lets skipstart until instructed otherwise
      echo "reverse" > skipstart
   fi

#  END_of_B___END_of_B___END_of_B___END_of_B___END_of_B___END_of_B___END_of_B___END_of_B___

#  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  propagation loop
   while (true)
   do
#increment runpointnumber
      cp runpointnumber $PROG_SCRDIR/temp533
      awk 'BEGIN {getline;i=$1+1;print i}' $PROG_SCRDIR/temp533 > runpointnumber
      rm $PROG_SCRDIR/temp533
      rm -f $PROG_SCRDIR/goingwell
      cp $TEMPORARY_DIR/g09.com $PROG_SCRDIR/g09.com
	 cp $PROG_SCRDIR/g09.com $COM_LOG_FILES/g09"$RUN_NUMBER".com
	 $G09_ROOT/bin/rung09 $PROG_SCRDIR/g09.com > $PROG_SCRDIR/g09.log
	 cp $PROG_SCRDIR/g09.log $COM_LOG_FILES/g09"$RUN_NUMBER".log
	 ((RUN_NUMBER++))
      grep 'Normal termination' $PROG_SCRDIR/g09.log > $PROG_SCRDIR/goingwell
      if (test -s $PROG_SCRDIR/goingwell) then
         awk -f $TEMPORARY_DIR/proganal $PROG_SCRDIR/g09.log >> $TEMPORARY_DIR/dynfollowfile
         mv olddynrun olderdynrun
         awk '/Input orientation/,/Distance matrix/ {print}' $PROG_SCRDIR/g09.log | awk '/   0   / {print}' > old
         cp $PROG_SCRDIR/g09.log olddynrun
         awk '/Input orientation/,/Distance matrix/ {print}' olderdynrun | awk '/   0   / {print}' > older
         awk -f $TEMPORARY_DIR/progdynb $PROG_SCRDIR/g09.log > g09.com
         rm -f old older 
      else
         cp $PROG_SCRDIR/g09.log $TEMPORARY_DIR/g09.log
         break
      fi
# kludge to do a side calculation of NMR using progcfour.  If ZMAT is there then it gets ran and renamed.
# creation of ZMAT is under the control of progdynb, which is controlled by keyword NMRcc in progdyn.conf
# decisions to be made: erase ZMAT at beginning?  what to do if cfour calc dies?
      if (test  -f ZMAT) then
         cp ZMAT $PROG_SCRDIR
         $PROG_SCRDIR/progcfour $TEMPORARY_DIR $PROG_SCRDIR
         
         mv ZMAT temp.ZMAT
         echo "generic one two three" `cat runpointnumber` "runisomer" `cat isomernumber` >> NMRlistcc
         awk '/Nuclear Magnetic Resonance/,/HF-SCF/ {if ($2=="C") print $1,$2,"Isotropic =",$3; if ($2=="H") print $1,$2,"Isotropic =",$3}' x.log >> NMRlistcc
      fi

# here is a cool link that lets you interupt the dynamics with a short job, then
# it automatically goes back to the dynamics  just make the file 'detour' and it
# will delete detour, run run.com, then go back to dynamics
      if (test  -f detour) then
         rm detour
         date >> $LOG_FILE
         cat run.com >> $LOG_FILE
         cp run.log temp.log
         $G09_ROOT/bin/rung09 $TEMPORARY_DIR/run.com > $TEMPORARY_DIR/run.log
         
      fi

#stop it all nicely by creating a nogo file
      if (test  -f nogo) then
         break
      fi

#figure out if this isomer is done - change in april 2013 is to move proganal call up from here
      rm -f $PROG_SCRDIR/tempdone
      tail -2 dynfollowfile | awk '/XXXX/ {print}' > $PROG_SCRDIR/tempdone
      if (test  -s $PROG_SCRDIR/tempdone) then
         if [ `awk '/reversetraj/ {if ($1=="reversetraj") print $2}' progdyn.conf` = "true" ]; then
            if [ `cat skipstart` = "reverse" ]; then
               rm -f skipstart
               rm -f geoPlusVel
               rm -f olddynrun
               rm -f olderdynrun
               a=`awk '{print $1}' isomernumber`
               mv traj traj$a
               mv dyn dyn$a
            fi
            if [ `cat skipstart` = "forward" ]; then
               echo reverserestart > skipstart
            fi
         else
            rm -f skipstart
            rm -f geoPlusVel
            rm -f olddynrun
            rm -f olderdynrun
            a=`awk '{print $1}' isomernumber`
            mv traj traj$a
            mv dyn dyn$a
         fi
         break
      fi
   done
#  END_of_C_Loop____END_of_C_Loop____END_of_C_Loop____END_of_C_Loop____END_of_C_Loop____END_of_C_Loop____

# We've got to break a second time to get out of this loop
# if we really want to quit.  Otherwise, it will start over
# at the top
   if (test  -f nogo) then
      break
   fi
   if (test  -s $PROG_SCRDIR/goingwell) then
      echo "starting a new point or a new direction"
   else
      break
   fi
done

   if (test  -f nogo) then
      break
   fi
   if (test  -s $PROG_SCRDIR/goingwell) then
      echo "starting a new point or a new direction2"
   else
      break
   fi
done
#exit 0
