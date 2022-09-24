prefix="$1"
grep 'Total energy' $prefix.log |cut -f2 -d= > ${prefix}_tot_energy.time
grep 'Kinetic energy' $prefix.log |cut -f2 -d= > ${prefix}_kinet_energy.time
grep 'Thermal energy' $prefix.log |cut -f2 -d= > ${prefix}_therm_energy.time
grep 'Gravitational energy' $prefix.log |cut -f2 -d= > ${prefix}_grav_energy.time
