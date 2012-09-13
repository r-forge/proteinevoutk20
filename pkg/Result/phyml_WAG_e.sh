nohup phyml -i rokasAA -d aa -q -m WAG -f e -v 0 -c 1 -u gtr_newick -o r --run_id WAG+e > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m WAG -f e -v e -c 1 -u gtr_newick -o r --run_id WAG+e+I > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m WAG -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id WAG+e+G > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m WAG -f e -v e -c 4 -a e -u gtr_newick -o r --run_id WAG+e+I+G > /dev/null & 
