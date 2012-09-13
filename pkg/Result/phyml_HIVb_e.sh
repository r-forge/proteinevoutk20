nohup phyml -i rokasAA -d aa -q -m HIVb -f e -v 0 -c 1 -u gtr_newick -o r --run_id HIVb+e > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVb -f e -v e -c 1 -u gtr_newick -o r --run_id HIVb+e+I > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVb -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id HIVb+e+G > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVb -f e -v e -c 4 -a e -u gtr_newick -o r --run_id HIVb+e+I+G > /dev/null & 
