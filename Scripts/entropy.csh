#!/bin/csh
set mols = "host G1 G2 G3 G3_Iso G4 G4_Iso G5 G6 G7"

set mols = "host"

foreach  mol ($mols)
         cd $mol
           #python ../../../bin/mcc_api_UA.py guest/guest.trr guest/guest.tpr guest/guest.log -1 -1 1 3 1 0.5 300 
           #python ../../../bin/mcc_api_UA.py complex/guest.trr complex/guest.tpr complex/guest.log -1 -1 1 3 1 0.5 300
           python ../../../bin/mcc_api_UA.py host.trr host.tpr host.log -1 -1 1 3 1 0.5 300
	   #python ../../../bin/mcc_api_UA.py water.trr water.tpr water.log -1 -1 1 3 1 0.5 300
           #cat Svib_M.txt Srot_UA.txt Strans_UA.txt | awk '$1!="END"' > Entropy.txt
        cd ..
end
