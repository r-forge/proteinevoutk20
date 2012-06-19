#!/bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -q medium*
#$ -M chaij@umail.iu.edu
#$ -m beas
#$ -S /bin/bash
#$ -pe threads 12
#$ -N B3_3
/data/apps/R/2.14.0/bin/R CMD BATCH bg3.3.R