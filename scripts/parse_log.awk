#!/bin/awk
#
#  awk -f murat/scripts/parse_g4_log.awk pet_001.log
#
BEGIN {
    file = fn;
    rc   = -1;
    n    = 1;
}

/SeedFinder::produce  event number/ {
    event = $5;
#   print $5;
}

/missed delta: sim.id/ {
    printf "event: %5i %s\n", event,$0;
}

END {
#    if (rc < 0) {
#	nevents = -1
#	cpu     = -1
#	real    = -1
#    }
#
#    printf "rc: %5i nevents: %6i CPUTime: %12.6f RealTime: %12.6f file: %s", rc, "nevents ", nevents,"CPUTime ", cpu, "RealTime: ", real, fn
#    printf "%5i %6i %12.6f %12.6f %s\n", rc,nevents,cpu,real,fn
}
