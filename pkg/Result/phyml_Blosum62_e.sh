phyml -i rokasAA -d aa -q -m Blosum62 -f e -v 0 -c 1 -u gtr_newick -o r --run_id Blosum62+e > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f e -v e -c 1 -u gtr_newick -o r --run_id Blosum62+e+I > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f e -v 0 -c 4 -a e -u gtr_newick -o r --run_id Blosum62+e+G > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f e -v e -c 4 -a e -u gtr_newick -o r --run_id Blosum62+e+I+G > /dev/null & 
