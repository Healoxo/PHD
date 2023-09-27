#!/bin/bash
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/spectrum.txt
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/Spectum_Script/spectrum.mdl
rm /home/thomas/Documents/0_Stage/My_pimms_prgm/all_spectrum.txt

Nh=0 #(0 18 19 20 21 22)
filter="Medium"
band="0.2-2"
step=0.01
T=6
python3 /home/thomas/Documents/0_Stage/My_pimms_prgm/Spectum_Script/spectrum_generator.py $Nh $fliter $band $step $T           #Genere 1 header et n ligne (model + output)
/home/thomas/Documents/0_Stage/PIMMS/pimms < /home/thomas/Documents/0_Stage/My_pimms_prgm/spectrum.txt > /home/thomas/Documents/0_Stage/My_pimms_prgm/output_spectum.txt 
python3 /home/thomas/Documents/0_Stage/My_pimms_prgm/Spectum_Script/plot_spectrum.py $band $step $Nh #Lui lis le spectrum.mdl
echo Fini
