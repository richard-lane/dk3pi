OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=foo.exe
fi

g++ src/bin_generated_decays.cpp -Ofast -o $OUTNAME `root-config --cflags --glibs` -lboost_filesystem

