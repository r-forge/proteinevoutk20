nohup phyml -i rokasAA -d aa -q -m HIVw -f e -v 0 -c 1 -u gtr_newick -o r --run_id HIVw+e > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVw -f e -v e -c 1 -u gtr_newick -o r --run_id HIVw+e+I > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVw -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id HIVw+e+G > /dev/null & 
 nohup phyml -i rokasAA -d aa -q -m HIVw -f e -v e -c 4 -a e -u gtr_newick -o r --run_id HIVw+e+I+G > /dev/null & 
