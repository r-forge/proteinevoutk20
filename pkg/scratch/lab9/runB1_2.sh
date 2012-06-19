#!/bin/bash
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null
#$ -q medium*
#$ -M chaij@umail.iu.edu
#$ -m beas
#$ -S /bin/bash
#$ -pe threads 12
#$ -N B1_2
/data/apps/R/2.14.0/bin/R CMD BATCH bg1.2.R