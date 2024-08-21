#!/usr/bin/env bash

inputfile=$1
output_directory=$2
filename=filelist
numfile=`wc -w $filename | awk '{print $1}'`

#To extract the filfilepath from inputfile
FilFilePath=`cat $inputfile | awk '{if (NR==39) print $2}'`
echo "$FilFilePath"

#To extract the fildir-path
fildir=`dirname "$FilFilePath"`

echo " $inputfile $output_directory $FilFilePath $fildir"

i=1
while [ $i -le $numfile ]
do
        scan=`cat $filename | head -$i | tail -1`
        echo "$scan of $numfile with ${scan} is getting analysed"


# Editing inpulfile for this scan
filfile=`echo ${scan}_IA_search.fil`
echo "$filfile"

newfilfile=`echo ${fildir}/${filfile}`
echo "$newfilfile"

cat $inputfile | head -28 >  temp.txt
echo "file $newfilfile" >> temp.txt
cat temp.txt >  $inputfile
rm -rf temp.pxt

currentdir=`pwd`

# Running Astro-Accelerate
echo "./astro-accelerate.sh $inputfile $output_directory"
./astro-accelerate.sh $inputfile $output_directory

DIR=$(dirname "$1")

cd $output_directory/output/GMRT.SPS.10m_newAA/
mkdir auto_output

# Running python sifting code
echo "python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat"
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 0.0 --ut 60.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 60.0 --ut 120.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 120.0 --ut 180.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 180.0 --ut 240.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 240.0 --ut 300.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 300.0 --ut 360.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 360.0 --ut 420.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 420.0 --ut 480.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 480.0 --ut 540.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 540.0 --ut 600.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 600.0 --ut 660.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 660.0 --SNR_threshold 6 #Change the threshold here

cd ${DIR}

# Running FETCH
echo "python auto_mod.py --pin $output_directory/output/GMRT.SPS.10m_newAA/FRB_detection_list.csv --pout $output_directory/output/GMRT.SPS.10m_newAA/  --fpath ${newfilfile}"
python auto_mod.py --pin $output_directory/output/GMRT.SPS.10m_newAA/FRB_detection_list.csv --pout $output_directory/output/GMRT.SPS.10m_newAA/  --fpath ${newfilfile}

rm -rf $output_directory/output/GMRT.SPS.10m_newAA/*fil

mv $output_directory/output/GMRT.SPS.10m_newAA $output_directory/output/GMRT.SPS.10m_newAA.$scan

cd $currentdir

i=$(( $i + 1 ))
done
