#!/bin/awk
#------------------------------------------------------------------------------
# debugging tool: print straw hits of a given MC particle
# example:
#         cat delta_finder_test1_event_0001.log | awk -v sim_id=35636 -v station=4 -f print_delta_seeds.awk
#--------------------------------------------------------------------------------------

BEGIN {
    doit = 0
}

/seed  / {
    st = $2
    if (st == station) station_ok = 1
    print $0
}

{
    if ($11 == sim_id) printf "%2i:%i %s\n",station,face,$0
}
