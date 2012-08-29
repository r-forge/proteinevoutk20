phyml -i rokasAA -d aa -q -m WAG -f m -v 0 -c 1 -u gtr_newick -o r --run_id WAG+m > /dev/null & 
 phyml -i rokasAA -d aa -q -m WAG -f m -v e -c 1 -u gtr_newick -o r --run_id WAG+m+I > /dev/null & 
 phyml -i rokasAA -d aa -q -m WAG -f m -v 0 -c 4 -a e -u gtr_newick -o r --run_id WAG+m+G > /dev/null & 
 phyml -i rokasAA -d aa -q -m WAG -f m -v e -c 4 -a e -u gtr_newick -o r --run_id WAG+m+I+G > /dev/null & 
