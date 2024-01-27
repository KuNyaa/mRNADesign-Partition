#!/bin/bash

# Go to ./src and run make
cd ./src
make
cd ..

# Go to ./tools/LinearPartition and run make
cd ./tools/LinearPartition
make
cd ../..

# Go to ./tools/LinearDesign and run make
cd ./tools/LinearDesign
make
cd ../..