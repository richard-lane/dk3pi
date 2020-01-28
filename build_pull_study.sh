OUTNAME=$1
if [[ -z $OUTNAME ]] ; then
    OUTNAME=bar.exe
fi

# Build with debug symbols without optimisation
g++ src/pull_study.cpp -Ofast -o $OUTNAME `root-config --cflags --glibs` \
-lboost_filesystem -Wall -Wextra -Wformat-security -Werror -O0 -g
