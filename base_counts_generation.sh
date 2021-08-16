#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA -b parameterB"
   echo -e "\t-a The path to the chromosome file directory."
   echo -e "\t-b The suffix for the 'base_counts_suffix.txt' file"
   exit 1 # Exit script after printing help
}

while getopts "a:b:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterB" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$parameterA"
echo "$parameterB"

touch ${parameterA}/base_counts_${parameterB}.txt

for f in ${parameterA}/*.fa; do

  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $f | \
   awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2,length($2)}' | \
   awk '{OFS="\t"}{print $1,$3}' >> ${parameterA}/base_counts_${parameterB}.txt

done

Rscript base_counts.R \
  ${parameterA} base_counts_${parameterB}.txt base_counts_${parameterB}.txt


#awk '{printf "%s\n", $2>$1".seq"}' file
