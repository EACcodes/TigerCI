#! /bin/python
#
# Checks the error for my Helium2 excited state test
# Necessary, because the weird output TigerCI produces
# for SDCI excited states is annoying to parse in bash

import sys

bench_state_1 = float(sys.argv[1])
bench_state_2 = float(sys.argv[2])
threshold     = float(sys.argv[3])
name          = sys.argv[4]

print name
energies = []
text = [ line for line in open("workdir/%s.output"%name, "r") ]
for line in text:
    if line.find('total energy') > 0 :
        tmp = line.strip().split()
        energies.append( float(tmp[-1]) )
        
diff_state_1 = abs(energies[0] - bench_state_1 )
diff_state_2 = abs(energies[1] - bench_state_2 )

print "final energies = %s %s " % ( energies[0], energies[1] )
print "benchmark      = %s %s " % ( bench_state_1, bench_state_2)
print "differences    = %s %s " % ( diff_state_1, diff_state_2 )

if ( abs(diff_state_1) < threshold and abs(diff_state_2) < threshold):
    print "Test is GOOD"
    sys.exit(0)
else:
    print "Test failed"
    f = open("workdir/FAILED","w")
    f.write("failed")
    f.close()
    sys.exit(1)

