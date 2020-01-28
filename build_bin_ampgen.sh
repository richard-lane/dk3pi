OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=foo.exe
fi

# Build with debug symbols without optimisation
g++ src/bin_generated_decays.cpp -o $OUTNAME `root-config --cflags --glibs` \
-lboost_filesystem -Wall -Wextra -Wformat-security -Werror -O0 -g
