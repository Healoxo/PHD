import matplotlib.pyplot as plt
import sys
Nh=0               #float(sys.argv[1])
filterr="Medium"   #str(sys.argv[2])
band="0.2-2"       #str(sys.argv[3])
lo,hi=band.split('-')
step=0.01          #str(sys.argv[4])
T=6                #float(sys.argv[5])
f=open("/home/thomas/Documents/0_Stage/My_pimms_prgm/spectrum.txt","w")
f.write(f"from xmm pn {filterr}  {band}\ninstrument flux ergs {band} unabsorbed")
f.write(f"\nmodel plasma {T} logt 1  {Nh}")
f.write(f"\noutput spectrum.mdl {lo} {hi} {step}")
f.close()
