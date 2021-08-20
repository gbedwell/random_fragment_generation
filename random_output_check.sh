#!/bin/bash

for file in random_fragments*.fa
do
  echo "$file $(grep -v '>' $file | grep -v '^[atgcATGCnN]*$' | wc -l) non-ATGCN line(s) $(grep -v '>' $file | grep '^$' | wc -l) empty line(s)"
done
