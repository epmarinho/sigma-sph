#!/usr/bin/bash
prefix="$1"
grep "Total energy" $prefix.log|cut -f4 -d' '  > Total_energy.his
grep Kinetic $prefix.log|cut -f4 -d' '  > Kinetic.his
grep Gravitational $prefix.log |cut -f4 -d' ' > Gravitational.his
