# dk3pi
Phase space binning and analysis for D -> K3pi events

NB might need to build with `BOOST_ROOT=<boost installation> cmake .. && make` if boost installation not found
NB; might need to build with
`cmake -DCMAKE_INCLUDE_PATH=/users/mh19137/lib/boost_1_73_0/ -DCMAKE_LIBRARY_PATH=/users/mh19137/lib/boost_1_73_0/stage/lib/ ..`
if on e.g. sc01 or a machine with the right version of boost installed outside of the usual place.

Needs Boost version > 1.66 (i think) for numerical integration

