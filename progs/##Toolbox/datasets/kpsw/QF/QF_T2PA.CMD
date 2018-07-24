/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

                             QFPV
                          Version 1.2

   This program is a revised version of QFPV, which is based on 'Testing
 for Common Trends,' (September, 1987) by J.H. Stock and M.W. Watson.
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

 Following is a list and descriptions of parameters defined in the
 program. Please read this and supply the appropriate values in the
 user definition area. (If the value is charater string, it should be
 enclosed in double quotes. Each value should also satisfy the condition
 given next to each description.)

 IFNAME : Name of the dataset
    ORG : Data organization [ (V/O) : by Variables or by Observations ]
  ALTIM : Total number of time periods
     BG : Starting period [ 1 <= BG <= ALTIM ]
     ED : Ending period [ BG < ED <= ALTIM ]
     NS : Number of data series
      K : Number of unit roots under the null [ K <= 6 ]
   MINL : Minimum AR lag [ MINL >= 1 ]
   MAXL : Maximum AR lag
   INCL : AR lag increment from MINL to MAXL
   OPTN : Option number of preprocessing
          (Please see the description of three options below.)
 OFNAME : Name of the output file

 ---------------------------------------------------------------------
    There are 3 preprocessing options. Please choose the number.
    1. Detrending the data.( i.e., regression on a constant and time)
    2. Removing means from the data. (i.e., regression on a constant)
    3. Do nothing. (i.e., raw data)
 ---------------------------------------------------------------------
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

@-------------------------------------------------------------------------@
@====================== User Definition Area =============================@
@-------------------------------------------------------------------------@

   ifname="\\kpsw\\data\\mpyr.prn";  ofname="\\kpsw\\qf\\qf_t2pa.out";
                  org="O";                    altim=168;
                    bg=24;                       ed=168;
                     ns=3;                          k=3;
                   minl=1;                       maxl=5;
                   incl=1;                       optn=1;

@ !!!!!!!! PLEASE MAKE SURE THAT EIG PROCEDURE HAS BEEN LOADED !!!!!!!!!! @

@ ::: Following command runs the program ::: @
#include \kpsw\qf\qfpv.prg;
