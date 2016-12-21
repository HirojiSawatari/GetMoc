#!/bin/csh -fv
#QSUB -l mpp_p=2 -q ded_1 -lt 40000 -lT 43000 -lM 60Mw -lm 60Mw -eo
# modified from mocdec1.csh. The difference is that here the data files
# are for present day configuration, moc is calc. for each basin
# Three files will be produced: eulerian, eddy-induced and total transport

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
  # set CASE = sb_570MaR5H2km_ni6
  # set CASE = sb_570H2kmni2_br1a
  # set CASE = sb_570MaR5PreOcn_ni
  # set CASE = sb_570H0kmni3br2_ct1
  # set CASE = sb_720H0kmni3_br1
  # set CASE = sb_720Mar0R5H2km_ice1
  # set CASE = sb_720H2kmice1_hy5
  # set CASE = 1870ctl4
   set CASE = test5
   
  #set DIR = $PROJECT/data_tmp_moc
   set DIR = /home/ygliu/archive/midHolocene/moc_tmp
  #set MSSIN = $SCRATCH/archive/${CASE}/ocn
  #set MSSIN = $PROJECT/archive/${CASE}/ocn
  #set MSSIN = $PROJECT/data_tmp/720Ma 
  #set MSSIN = $PROJECT/data_tmp/1870atm_test 
   set MSSIN = /home/ygliu/archive/midHolocene

  #set MSSOUT = $PROJECT/data_tmp 
  set MSSOUT = /home/ygliu/archive/midHolocene/moc

  set DATADIR = ${MSSOUT}

  set yr1 = 1920
  set yr2 = 1969

  #don't add .nc to the file name
  ;set FILE0 = $CASE.meanOcn.yrs${yr1}-${yr2}ISOP
  set FILE0 = $CASE.meanOcn.yrs${yr1}-${yr2}
  set FILE  = ${FILE0}mean.nc

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
   cp $HOME/analys_soft/ocn/moc/GNUmakefile_basin GNUmakefile 
   cp $HOME/analys_soft/ocn/moc/*.f90 .

  gmake

  # produce the annual mean 
  rsync -avz $MSSIN/${FILE0}.nc .
  ncra ${FILE0}.nc ${FILE}

  #2)Run Program:
  #make input file moc.in:
  #change nphi if have different grid. It should be smaller than the number
  # of latitude used in the pop model. In my test, nphi=71 works for gx3v5
  #resolution, but 76 doesn't
  cat > moc.in << EOFINPUT1
 &io_info
  infile = '$FILE'
/
 &grid_info
  nphi = 71 
  phi_start = -90.0
  phi_end =  90.0
/
EOFINPUT1
               
  #Run moc:
  ./basins_pop_moc_eul < moc.in
  ./basins_pop_moc_ei < moc.in

  
  #calculate the total
  #ncbo -op_typ=add TVTS_${FILE} TVTSei_${FILE} -o TVTStot_${FILE}
  ncflint -w 1,1 TVTS_${FILE} TVTSei_${FILE} TVTStot_${FILE}
  rsync -avz  TVTS_${FILE}  $MSSOUT/${CASE}.TVTS.${yr1}_cat_${yr2}.nc
  rsync -avz  TVTSei_${FILE}  $MSSOUT/${CASE}.TVTSei.${yr1}_cat_${yr2}.nc
  rsync -avz  TVTStot_${FILE}  $MSSOUT/${CASE}.TVTStot.${yr1}_cat_${yr2}.nc
  rm TVTS*.nc
  rm ${FILE0}.nc ${FILE}
        
exit(0)
