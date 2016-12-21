#!/bin/csh -fv
#QSUB -l mpp_p=2 -q ded_1 -lt 40000 -lT 43000 -lM 60Mw -lm 60Mw -eo

set echo
set verbose
#################################################################
# This script makes the barotropic streamfunction and
# meridional overturning data using fortran code.
# Writes netcdf files.
# Loop over Centuries, decades, years.
#################################################################
#  Script by E. Brady
#  pop_moc modified from Hu, Strand, Tony
#  Merid. Overturning idl modified from G. Danabasoglu 
#################################################################
# script modified by c shields... no control case, only compute data
# plotting not added yet 
#################################################################

# DIR is working directory:
  # set CASE = sb_570MaR5H0km_ni3
  # set CASE = sb_570H2kmni2_br1a
  # set CASE = sb_570MaR5PreOcn_ni
  set CASE = sb_570H0kmni3_ct1
  set DIR = $PROJECT/data_tmp_moc
  set MSSIN = $SCRATCH/archive/${CASE}/ocn
  #set MSSIN = $PROJECT/archive/${CASE}/ocn
  set MSSOUT = $PROJECT/data_tmp 
  set DATADIR = ${MSSOUT}

#To loop over years, decades, centuries...
set CENTURIES = ( 40 )                  # centuries to loop over
set DECADES   = ( 0  ) # decades to loop over
# set YEARS     = ( 1 ) # years to loop over in 1yr increments
 set YEARS     = ( 0 1 2 3 4 5 6 7 8 9) # years to loop over in 1yr increments
#######################################################################
#Don't Need to Change below here if following csm conventions
#    for file names and MASSSTORE data locations:
------------------------------------------------------
# cd to a temporary directory; get needed files etc.
  if !(-d $DIR)  mkdir -p $DIR
  cd $DIR 
#-------------------------------------------------------------
#1)Compute Meridional Overturning Streamfunction (basins)
#Compile program:
 cp $HOME/snowball_ccsm3/analys_soft/ocn/moc/GNUmakefile .
 cp $HOME/snowball_ccsm3/analys_soft/ocn/moc/*.f90 .

 gmake

   set endyr = $YEARS[$#YEARS]
   set begyr = $YEARS[1]   
   set enddec = $DECADES[$#DECADES]
   set begdec = $DECADES[1]     
   set endcen = $CENTURIES[$#CENTURIES]
   set begcen = $CENTURIES[1]     # single year name  (for ann output filename)
##############
#  because in CCSM: 1st year of first decade of first century = 1:
   if  ( $begcen == 00 ) then
     if ($begdec == 0 && $begyr == 0 ) set begyr = $YEARS[2]
   endif
##############
   setenv Decfinalb `printf "%02d%01d%01d\n" $begcen $begdec $begyr`
   setenv Decfinale `printf "%02d%01d%01d\n" $endcen $enddec $endyr`
   setenv Decfinal  ${Decfinalb}_cat_${Decfinale}

   echo $Decfinal
   if (-f $CASE.pop.h.TVTS.$Decfinal.nc ) goto 1002 

   foreach CC ( $CENTURIES )    # loop over centuries
       foreach DD ( $DECADES )     # loop over decades
           set begyr = $YEARS[1]
           # if ($CC == 00 && $DD == 0 && $#YEARS > 1) set begyr = $YEARS[2]
           echo $#YEARS
           set endyr = $YEARS[$#YEARS]
        
           setenv FNfinalb `printf "%02d%01d%01d\n" $CC $DD $begyr`
           setenv FNfinale `printf "%02d%01d%01d\n" $CC $DD $endyr`
           setenv FNfinal  ${FNfinalb}_cat_${FNfinale}
        
           if (-f $CASE.pop.h.TVTS.$FNfinal.nc ) goto 1001 
         
           foreach YY ( $YEARS )      # loop over years
         
               setenv FN `printf "%02d%01d%01d\n" $CC $DD $YY`
               echo starting $FN `date`
               set FILE = $CASE.pop.h.$FN.nc 
               
               #####
               # added by ygliu, August 15, 2011
               # To copy and produce the annual mean for year FN
               set FILES = $CASE.pop.h.$FN
               ncra $MSSIN/${FILES}-??.nc $MSSIN/${FILE}
               rsync -avz $MSSIN/${FILE} .
               #####
             
               # rsync -avz $MSSIN/$FILE .
               if !(-f $FILE ) goto 1000 
             
               #2)Run Program:
               #make input file moc.in:
               #change nphi if have different grid
               cat > moc.in << EOFINPUT1
 &io_info
  infile = '$FILE'
/
 &grid_info
  nphi = 181 
  phi_start = -90.0
  phi_end =  90.0
/
EOFINPUT1
               
               #Run moc:
               pop_moc < moc.in
               pop_moc_ei < moc.in
               
               ###############################
               # skipping over year 0 decade 0
               1000:
               ###############################
           end #Years
           # concat into decade bundles and archive:
           ncrcat TVTS_$CASE.pop.h.*.nc $CASE.pop.h.TVTS.$FNfinal.nc 
           ncrcat TVTSei_$CASE.pop.h.*.nc $CASE.pop.h.TVTSei.$FNfinal.nc 
           rsync -avz  $CASE.pop.h.TVTS.$FNfinal.nc  $MSSOUT/$CASE.pop.h.TVTS.$FNfinal.nc
           rsync -avz  $CASE.pop.h.TVTSei.$FNfinal.nc  $MSSOUT/$CASE.pop.h.TVTSei.$FNfinal.nc
           rm TVTS*.nc
           rm $CASE.pop.h.0*.nc
        
           ######################################
           # skipping over existing decadal files 
           1001:
           ######################################
       end #Decades
   end #Centuries
  
   if ($Decfinal !~ $FNfinal) ncrcat $CASE.pop.h.TVTS.*.nc $CASE.pop.h.TVTS.$Decfinal.nc 
   if ($Decfinal !~ $FNfinal) ncrcat $CASE.pop.h.TVTSei.*.nc $CASE.pop.h.TVTSei.$Decfinal.nc 
  
   ########################################
   # skipping over creating long timeseries 
   1002:
   ########################################
  
   mv $CASE.pop.h.TVTS.$Decfinal.nc $DATADIR/$CASE.TVTS.$Decfinal.nc
   mv $CASE.pop.h.TVTSei.$Decfinal.nc $DATADIR/$CASE.TVTSei.$Decfinal.nc
   rm *.nc
  
exit(0)
