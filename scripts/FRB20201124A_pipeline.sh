#!/usr/bin/env bash

inputfile=$1
output_directory=$2
filename=filelist
numfile=`wc -w $filename | awk '{print $1}'`

#To extract the filfilepath from inputfile
FilFilePath=`cat $inputfile | awk '{if (NR==33) print $2}'`
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
#newfilfile=`echo ${FilFilePath}`

echo "$newfilfile"

cat $inputfile | head -32 >  temp.txt
echo "file $newfilfile" >> temp.txt
cat temp.txt >  $inputfile
rm -rf temp.pxt

currentdir=`pwd`

# Running Astro-Accelerate
echo "./astro-accelerate.sh $inputfile $output_directory"
./astro-accelerate.sh $inputfile $output_directory

DIR=$(dirname "$1")

cd $output_directory/output/FRB20201124A/
mkdir auto_output

# Running python sifting code
echo "python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat"
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 0.0 --ut 300.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 300.0 --ut 600.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 600.0 --ut 1200.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 1200.0 --ut 1800.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 1800.0 --ut 2400.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 2400.0 --ut 3000.0 --SNR_threshold 6 #Change the threshold here
python ${DIR}/spsplotii_frb_detect_py3.py global_peaks.dat --lt 3000.0 --ut 3600.0 --SNR_threshold 6 #Change the threshold here

cd ${DIR}

# Running FETCH
echo "python auto_mod.py --pin $output_directory/output/FRB180301/FRB_detection_list.csv --pout $output_directory/output/FRB180301/auto_output/ --fpath ${newfilfile}"
python auto_mod.py --pin $output_directory/output/FRB20201124A/FRB_detection_list.csv --pout $output_directory/output/FRB20201124A/auto_output/ --fpath ${newfilfile}

rm -rf $output_directory/output/FRB20201124A/auto_output/*fil

mv $output_directory/output/FRB20201124A $output_directory/output/$scan

cd $currentdir

i=$(( $i + 1 ))
done
