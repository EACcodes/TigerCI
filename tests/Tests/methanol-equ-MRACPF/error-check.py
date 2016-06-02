import sys

# inputs
bench_ref      = float(sys.argv[1])
bench_total    = float(sys.argv[2])
output         = sys.argv[3]
tol            = float(sys.argv[4])

# read the input file 
text = [ line.strip() for line in open(output, 'r')]

energy = [ line for line in text if line.find('total energy (electronic + nuclear) is :') >= 0 ]
energy = (energy[0].split())[-1]
energy = float(energy)


refEnergy = [ line for line in text if line.find('the reference energy is (elec + nuclear)') >= 0]
refEnergy = (refEnergy[0].split())[-1]
refEnergy = float(refEnergy)

ref_check   = abs(refEnergy - bench_ref) < tol
total_check = abs(bench_total - energy)  < tol

print 'Calculated reference energy = ' , refEnergy
print 'Benchmark  reference energy = ' , bench_ref
print 'Difference  = ' , abs(refEnergy - bench_ref)
print ' ' 
print 'Calculated total energy = ', energy
print 'Benchmark total energy  = ' , bench_total
print 'Difference = ' , abs(bench_total - energy)

if ref_check and total_check:
    print "Test is GOOD"
    sys.exit(0)
else:
    print 'Test Failed' 
    sys.exit(1)
