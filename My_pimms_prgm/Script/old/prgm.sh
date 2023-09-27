#!/bin/bash
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/all.txt
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt
Nh_values=(0 18 19 20 21 22)
#read -p "Enter count rate: " cts
#read -p "Enter filter(Thin,Medium,Thick): " filter
#read -p "Choice of the band (0.2-2 or 2-12): " band
cts=1e-2
filter="Medium"
band="0.2-2"
for Nh in "${Nh_values[@]}"
do
    python3 input_generator.py $Nh $cts $filter $band
    /home/thomas/Documents/0_Stage/PIMMS/pimms < /home/thomas/Documents/0_Stage/My_pimms_prgm/input.txt > /home/thomas/Documents/0_Stage/My_pimms_prgm/output.txt 
    python3 output_convertissor.py 
done
python3 read_all_output.py $cts
