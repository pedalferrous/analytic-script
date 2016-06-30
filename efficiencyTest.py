from time import *
from numpy import *
from itertools import *
import ast

f = open("test_list.txt", "r")
test_list = ast.literal_eval(f.read())
f.close()

# print test_list[0][0]
# print len(test_list)

t0 = clock()
for x in xrange(5000):
	mean0 = (sum(test_list[j][k][l]
	    for j in range(len(test_list))
	    for k in range(len(test_list[0]))
	    for l in range(len(test_list[0][0])))
	    /(len(test_list)*len(test_list[0])*len(test_list[0][0])))
t1 = clock()

t2 = clock()
for x in xrange(5000):
	mean1 = mean(list(chain.from_iterable(chain.from_iterable(test_list))))
t3 = clock()

# mean0 = 0
# mean1 = 0
# t0 = 0
# t1 = 0
# t2 = 0
# t3 = 0

print "__________ RESULTS __________"
print "Found mean of " + str(mean0) + " in " + str((t1-t0)/5.0) + " ms"
print "Found mean of " + str(mean1) + " in " + str((t3-t2)/5.0) + " ms"
print ""
