#!/bin/bash

echo "This is the list of the all the ITHACAoutput folder in your ITHACA-FV directory:"
echo ""
find . -name "ITHACAoutput"
echo ""
read -p "Are you sure you want to delete all of them (y/n)?" choice
case "$choice" in 
  y|Y ) find . -name "ITHACAoutput" -exec rm -r "{}" \;;;
  n|N ) echo "You decided not to remove them.";;
  * ) echo "Invalid answer, Type y or n";;
esac
#find . -name "ITHACAoutput" -delete 
