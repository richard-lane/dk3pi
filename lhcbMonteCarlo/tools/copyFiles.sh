#!/bin/bash
# Copy root files created from a ganga job to a (new) directory called files/
# Then squish them all together with ROOT's hadd command thing

# Some stuff is hard coded like the number of file

set -e

mkdir files


for i in {0..22}
do
    cp ~/gangadir/workspace/rilane/LocalXML/17/$i/output/DVntuple.root files/DVntuple${i}.root
done

hadd bigRootFile.root files/*.root

