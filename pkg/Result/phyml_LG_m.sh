nohup phyml -i rokasAA -d aa -q -m LG -f m -v 0 -c 1 -u gtr_newick -o r --run_id LG+m > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m LG -f m -v e -c 1 -u gtr_newick -o r --run_id LG+m+I > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m LG -f m -v 0 -c 4 -a e -u gtr_newick -o r --run_id LG+m+G > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m LG -f m -v e -c 4 -a e -u gtr_newick -o r --run_id LG+m+I+G > /dev/null & 
