#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -a parameterA -b parameterB -c parameterC"
   echo -e "\t-a The path to the large multi-entry fasta file."
   echo -e "\t-b The filename for the output one-entry-per-line file."
   echo -e "\t-c The directory to store the individual chromosome fasta files."
   exit 1 # Exit script after printing help
}

while getopts "a:b:c:" opt
do
   case "$opt" in
      a ) parameterA="$OPTARG" ;;
      b ) parameterB="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$parameterA" ] || [ -z "$parameterB" ] || [ -z "$parameterC" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$parameterA"
echo "$parameterB"
echo "$parameterC"


awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ${parameterA} > ${parameterB}

mkdir ${parameterC}

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > ${parameterC}/${outfile}
    else
        echo $line >> ${parameterC}/${outfile}
    fi

done < ${parameterB}

rm $parameterB
