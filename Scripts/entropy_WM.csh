#!/bin/csh
set mols = "host G1 G2 G3 G3_Iso G4 G4_Iso G5 G6 G7"

set mols = "host"

foreach  mol ($mols)
         cd $mol
           #mcc_gromacs.py --tpr guest/guest.tpr --trr guest/guest.trr -o guest/guest_WM.log -b 1 -e 100000000 --temper 300 --molecule -v 5 --fscale 1 --tscale 1
	   #mcc_gromacs.py --tpr complex/guest.tpr --trr complex/guest.trr -o complex/guest_WM.log -b 1 -e 100000000 --temper 300 --molecule -v 5 --fscale 1 --tscale 1 
           #mcc_gromacs.py --tpr complex/host.tpr --trr complex/host.trr -o complex/host_WM.log -b 1 -e 100000000 --temper 300 --molecule -v 5 --fscale 1 --tscale 1
	  #mcc_gromacs.py --tpr complex/complex.tpr --trr complex/complex.trr -o complex/complex_WM.log -b 1 -e 100000000 --temper 300 --molecule -v 5 --fscale 1 --tscale 1
	  mcc_gromacs.py --tpr host.tpr --trr host.trr -o host_WM.log -b 1 -e 100000000 --temper 300 --molecule -v 5 --fscale 1 --tscale 1 
          #cat Svib_M.txt Srot_UA.txt Strans_UA.txt | awk '$1!="END"' > Entropy.txt
        cd ..
end
