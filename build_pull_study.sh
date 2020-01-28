OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=pull.exe
fi

# Build without debug symbols, including optimisation
g++ src/pull_study.cpp -o $OUTNAME `root-config --cflags --glibs` \
-lboost_filesystem -Wall -Wextra -Wformat-security -Werror -O3 -g
