phyml -i rokasAA -d aa -q -m Blosum62 -f m -v 0 -c 1 -u gtr_newick -o r --run_id Blosum62+m > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f m -v e -c 1 -u gtr_newick -o r --run_id Blosum62+m+I > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f m -v 0 -c 4 -a e -u gtr_newick -o r --run_id Blosum62+m+G > /dev/null & 
 phyml -i rokasAA -d aa -q -m Blosum62 -f m -v e -c 4 -a e -u gtr_newick -o r --run_id Blosum62+m+I+G > /dev/null & 
