import sys
import pickle

if len(sys.argv) != 3:
    print "Usage: {} <pkl> <dom>\n";
    exit()




d2ind = pickle.load(open(sys.argv[1],'r'))
ind2d = {i:d for d,i in  d2ind['ObjectA'].items()}

with open('{}_status'.format(sys.argv[2]),'r') as openfile:
    bestV = float('-inf')
    bestI = 0
    i = 0
    for line in openfile:
        dat = line.replace(':','').split(' ')

        dat = dat[4:]
        neg = False
        for d in dat:
            if d != "":
                n = float(d)
                if n < 0:
                    if neg:
                        v = n
                        break
                else:
                    neg = True
        if v > bestV:
            bestV = v
            bestI = i
        i += 1
print bestV,bestI
bestCategories = []
catOccurrences = {i:{} for i in d2ind['ObjectA']}
with open('{}_dom0'.format(sys.argv[2]),'r') as openfile:
    i = -1
    for line in openfile:
        if True: #i == bestI:
            cats = line.rstrip().split(' ')
            if i == bestI:
                bestCategories = cats
            categories = {}
            for ii in range(len(cats)):
                cat = cats[ii]
                if cat not in categories:
                    categories[cat] = []
                categories[cat].append(ind2d[ii])
            for cat in categories.values():
                for c1 in cat:
                    for c2 in cat:
                        if c1 != c2:
                            if c2 not in catOccurrences[c1]:
                                catOccurrences[c1][c2] = 0
                            catOccurrences[c1][c2] += 1
                            
            #for cat in categories.values():
            #    print cat


        i += 1
#print cats
#print d2ind['ObjectA']

ft = 0
fd = 0
max_ = float('-inf')
for cat in catOccurrences:
    total = 0
    for c in catOccurrences[cat]:
        if catOccurrences[cat][c] > max_:
            max_ =  catOccurrences[cat][c]
        total += catOccurrences[cat][c]
    fd += len(catOccurrences[cat])
    ft += total
for cat in catOccurrences:
    total = 0
    for c in catOccurrences[cat]:
        total += catOccurrences[cat][c]
    mean = float(ft)/float(fd)
    for c in catOccurrences[cat]:
        if catOccurrences[cat][c] > (mean + max_)*0.5:
            pass # print cat, c,catOccurrences[cat][c],(mean+ max_)*0.5
        

categories = {}
cats = bestCategories
print bestCategories
for i in range(len(cats)):
    cat = cats[i]
    if cat not in categories:
        categories[cat] = []
    categories[cat].append(ind2d[i])

for cat in categories.values():
    print cat
