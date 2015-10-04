#!/bin/bash
for file in $(echo $1 | tr "," " "); 
do
    echo $file
done


