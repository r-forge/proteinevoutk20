phyml -i rokasAA -d aa -q -m MtREV -f e -v 0 -c 1 -u gtr_newick -o r --run_id MtREV+e > /dev/null & 
 phyml -i rokasAA -d aa -q -m MtREV -f e -v e -c 1 -u gtr_newick -o r --run_id MtREV+e+I > /dev/null & 
 phyml -i rokasAA -d aa -q -m MtREV -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id MtREV+e+G > /dev/null & 
 phyml -i rokasAA -d aa -q -m MtREV -f e -v e -c 4 -a e -u gtr_newick -o r --run_id MtREV+e+I+G > /dev/null & 
