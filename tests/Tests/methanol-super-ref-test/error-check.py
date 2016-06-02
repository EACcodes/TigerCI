import sys

# inputs
bench_total    = float(sys.argv[1])
output         = sys.argv[2]
tol            = float(sys.argv[3])

# read the input file 
text = [ line.strip() for line in open(output, 'r')]

energy = [ line for line in text if line.find('total energy (electronic + nuclear) is :') >= 0 ]
energy = (energy[0].split())[-1]
energy = float(energy)

total_check = abs(bench_total - energy)  < tol

print 'Calculated total energy = ', energy
print 'Benchmark total energy  = ' , bench_total
print 'Difference = ' , abs(bench_total - energy)

if  total_check:
    print "Test is GOOD"
    sys.exit(0)
else:
    print 'Test Failed' 
    sys.exit(1)
