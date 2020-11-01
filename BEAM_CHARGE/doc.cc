10/05/2020
I  am going to calculate the beam charge for H/T data and T/T data eor the E12-17-003 experiment.
Beam charge is used to calculate the cross section
Original Root Files are required to calculate the beam charge.
The current cut can be used to select the beam current that the current is > 10 uA or any value can be selected.

The original code is given by Shujie
/work/halla/triton/shujie/get_charge.C


The original code is modified as BY DEVILAL ADHIKARI
code: run_charge.C
shell: find_charge.sh
list: run_list.list    ---> list of the runs to calculate the beam charge

to run the  code ./find_charge.sh

this code will produce a final output file as "outcharge.csv" which contains the following information.

run #,   beam time in min,beam current, beam charge in coulumb and total events in the run

  111221(run #),   27.895 ( beam time min), 22.276 (beam current),0.0372(beam cherge coulumb), 519772 (total events in the run)  



 The value of the beam current or beam charge should not be negative. If for some runs the value of the beam charge is negative, it must be corected. The way to correct is the following
  a. Add the total beam charge for the runs having +ve value of beam curent. Lets denote this value by T1
  1.  add the beam charge for about 25  runs or more  having the positive value of the beam current and calculate the beam charge per event or the beam charge per trigger. Lets say this quantity A
  2. Sum all of the events for the runs having the negative value of beam current and calculate the total events. Lets say this quantity B
  3. Now Multiply A*B wil give the beam charge for the runs having the negative value of current.Lets say this quantity C
  4. Add the C +T1   === total value of beam curents for all of the runs.
