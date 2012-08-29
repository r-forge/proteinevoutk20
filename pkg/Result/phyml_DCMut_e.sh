phyml -i rokasAA -d aa -q -m DCMut -f e -v 0 -c 1 -u gtr_newick -o r --run_id DCMut+e > /dev/null & 
 phyml -i rokasAA -d aa -q -m DCMut -f e -v e -c 1 -u gtr_newick -o r --run_id DCMut+e+I > /dev/null & 
 phyml -i rokasAA -d aa -q -m DCMut -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id DCMut+e+G > /dev/null & 
 phyml -i rokasAA -d aa -q -m DCMut -f e -v e -c 4 -a e -u gtr_newick -o r --run_id DCMut+e+I+G > /dev/null & 
