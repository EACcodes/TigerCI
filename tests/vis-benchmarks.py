#! /usr/bin/python
#
# visualize the verbose benchmark output
import sys
import numpy as np
import matplotlib.pyplot as plt
from pyparsing import Word, nums, alphanums, CaselessLiteral, Literal, Optional, Combine, OneOrMore 

"""

    Grammar Rules


"""


# floating point grammar 
point = Literal('.')
e = CaselessLiteral('E')
plusorminus = Literal('+') | Literal('-')
number = Word(nums+".") 
integer = Combine( Optional(plusorminus) + number )
floatnumber = Combine( integer +
                       Optional( point + Optional(number) ) +
                       Optional( e + integer )
                     )

# benchmark grammar
title        = Word(alphanums + "_" + "-")
stars       = OneOrMore(Literal("*"))
header      = title("name") + stars + Literal("testing the final energy")
result      = Literal("Test is GOOD")("pass") | Literal("Test failed")("fail")
single_root = ( Literal("final energy =") + Optional(floatnumber)("energy") + 
                Literal("benchmark    =") + floatnumber + 
                Literal("difference   =") + floatnumber("diff") ) 

multi_root  = ( Literal("Root") + number("root") + 
                Literal("Calculated total energy =") + Optional(floatnumber("energy")) + 
                Literal("Benchmark total energy  =") + floatnumber + 
                Literal("Difference =") + floatnumber("diff") )

entire_single_root = header + single_root + result
entire_multi_root  = header + OneOrMore(multi_root) + result

# timings grammar
header = Literal("INFO COMING FOR") + title("name")
data_points = Literal('Number of data points in "') + title + Literal('/bench.dat" =') + floatnumber('n')
mean = Literal("Arithmetic mean (average) =") + floatnumber('avg')
deviation = Literal("Standard Deviation =") + floatnumber('dev')        
timing_block = header + data_points + mean + deviation

"""

    Main program 
    
    
"""

def timings_parse(text):
    # parser grammar
    number = Word(nums + ".")
    test_name = Word(alphanums+"_"+"-")
    header = Literal("INFO COMING FOR") + test_name("name")
    data_points = Literal('Number of data points in "') + test_name + Literal('/bench.dat" =') + number('n')
    mean = Literal("Arithmetic mean (average) =") + number('avg')
    deviation = Literal("Standard Deviation =") + number('dev')        
    result_block = header + data_points + mean + deviation

    # read the data file
    data = {}
    matches = result_block.scanString(text)
    for match in matches:
        name =  match[0].name
        n = match[0].n
        avg = float(match[0].avg)
        dev = float(match[0].dev)
        data[name] = [n, avg, dev]

    # graph the data
    names = []
    values = []
    dev = []
    for point in data:
        names.append(point)
        values.append(data[point][1])
        dev.append(data[point][2])
    y = np.arange(len(names))
    
    plt.barh(y, values, xerr=dev, alpha=0.4)
    plt.yticks(y, names)
    plt.xlabel("Time (s)")
    plt.title("Benchmark Timings")
    plt.ylim((0,len(names)))
    (plt.gcf()).tight_layout(pad=1.09)
    plt.show()
    
def process_accuracies(matches):
    benchmarks = {}
    for match, start, stop in matches:
        name = match.name
        if not name in benchmarks: 
            benchmarks[name] = { 'passed':0, 'crashed':0, 'wrong':0 }
        if 'pass' in match:
            benchmarks[name]['passed'] += 1
        elif 'energy' in match:
            # if we didn't pass and the benchmarks found a total energy we == wrong answer
            # otherwise it just crashed (and the script won't find a total energy)
            benchmarks[name]['wrong'] += 1
        else:
            benchmarks[name]['crashed'] += 1
    return benchmarks
    
def accuracy_parse(text):
    # read and process the mail
    benchmarks = {}
    matches = entire_single_root.scanString(text)
    single_bench = process_accuracies(matches)
    matches = entire_multi_root.scanString(text)
    multi_bench =  process_accuracies(matches)

    # consolidate all the benchmarks for plotting
    names   = single_bench.keys() + multi_bench.keys()
    passed  = []
    crashed = []
    wrong   = []
    for name in names:
        if name in single_bench:
            target = single_bench[name]
        else:
            target = multi_bench[name]
        passed.append(target['passed'])
        crashed.append(target['crashed'])
        wrong.append(target['wrong'])


    # lets plot them
    y = np.arange(len(names))
    h = 0.2
    #plt.barh(y-0.3, passed, align='center', color='black', height=h)
    plt.barh(y-0.1, crashed, align='center', color='orange', height=h)
    plt.barh(y+0.1, wrong, align='center', color='red', height=h)
    plt.yticks(y, names)
    plt.xlabel("Number of Benchmarks")
    plt.title("Number of Failed/Crashed Benchmarks in Red/Orange")

    # fun manipulation of things to make it pretty ....
    ax = plt.gca()                                         # current axis
    fig = plt.gcf()                                        # current figure 

    x_max = max( max(crashed), max(wrong))
    plt.axis([0, x_max, 0, len(names)])                    # set the x axis to the minimum necessary size
    ax.set_xticks(np.arange(x_max)+1)                      # set integer values for the x axis ticks
    fig.tight_layout(pad=1.09)
    plt.show()
    
if __name__ == "__main__":
    try:
        file_name = sys.argv[1]
    except IndexError:
        print "This script requires the name of the text file containing the benchmarks output"
        sys.exit(0)
    try:
        text = open(file_name, "r").read()
    except Exception as e:
        print "Failure on opening/reading of the file %s" % file_name
        print e
    accuracy_parse(text)
    timings_parse(text)


