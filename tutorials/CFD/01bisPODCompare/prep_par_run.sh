#!/bin/bash


for i in {0..3}; do
  DEST="processor$i"
  echo "Copying files to $DEST"
  mkdir -p $DEST/system
  cp -r system/controlDict "$DEST/system/controlDict"
  cp -r system/ITHACAdict "$DEST/system/ITHACAdict"
  cp -r system/ITHACAPODdict "$DEST/system/ITHACAPODdict"
  cp -r constant/transportProperties "$DEST/constant/transportProperties"
done
