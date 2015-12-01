import sys
import pickle
import numpy as np
counts = {}
opposite = {'U':'D',
             'D':'U',
             'L':'R',
             'R':'L'}
kk = 0

effectCounts = {}
totalEffects = 0
jointProb = {}
filtered = {'Fireball'}
threshold = float(sys.argv[2])
objects = set()
with open(sys.argv[1],'r') as openfile:
    types = openfile.readline().rstrip().split(',')
    ind2type = {i:t for i,t in zip(range(len(types)),types)}
    typeSets = {t:set() for t in types}
    for line in openfile:
        dat = line.rstrip().split(',')
        oA = dat[0]
        oB = dat[1]
        A2B = dat[2]
        if oA not in counts:
            counts[oA] = {}
        if oB not in counts:
            counts[oB] = {}
        t = dat[3]
        objects.add(oA)
        objects.add(oB)
        objects.add(dat[4])
        
        ta = dat[5]
        dv = dat[6]
        if t != 'None' and not (oA in filtered or oB in filtered or dat[4] in filtered):
            kk += 1
            if dat[0] == dat[4]:
                s = 'A'
            if dat[1] == dat[4]:
                s = 'B'
            if oA not in counts:
                counts[oA] = {}
            acause = (oB,A2B)
            effect = (t,dat[4],ta,dv)
            if acause not in counts[oA]:
                counts[oA][acause] = {}
            
            if effect not in counts[oA][acause]:
                counts[oA][acause][effect] = 0
                
            counts[oA][acause][effect] += 1
            bcause = (oA,opposite[A2B])
            if oB not in counts:
                counts[oB] = {}
                
            if bcause not in counts[oB]:
                counts[oB][bcause] = {}
                
            if effect not in counts[oB][bcause]:
                counts[oB][bcause][effect] = 0
                
            counts[oB][bcause][effect] += 1

            totalEffects += 2
            if effect not in effectCounts:
                effectCounts[effect] = 0
            effectCounts[effect] += 2
            ceA = (acause,effect)
            ceB = (bcause,effect)
            if ceA not in jointProb:
                jointProb[ceA] = 0
            if ceB not in jointProb:
                jointProb[ceB] = 0
            jointProb[ceA] += 1
            jointProb[ceB] += 1
jointProb = {ce:float(c)/float(totalEffects) for ce,c in jointProb.items()}
effectProbabilities = {e:float(c)/float(totalEffects) for e,c in effectCounts.items()}
pmi = {}
npmi = {}
for o in counts:
    pmi[o] = {}
    npmi[o] = {}
    for cause in counts[o]:
        oEffects = 0
        for effect in counts[o][cause]:
            oEffects += counts[o][cause][effect]
        for effect in counts[o][cause]:
            condP = float(counts[o][cause][effect])/float(oEffects)
            pmi[o][(cause,effect)] = np.log(condP/effectProbabilities[effect])
            npmi[o][(cause,effect)] = np.log(condP/effectProbabilities[effect])/-np.log(jointProb[(cause,effect)])
p = []
npi = []

o = 'QBLOCK'
print o, pmi[o]
ceClusters = {}

for o in pmi:
    for ce in pmi[o]:
        p.append(pmi[o][ce])
        npi.append(npmi[o][ce])
        if pmi[o][ce] > threshold:
            if ce not in ceClusters:
                ceClusters[ce] = []
            ceClusters[ce].append(o)
        #print o,ce,pmi[o][ce],npmi[o][ce]
import matplotlib.pyplot as plt

clusters = set()
cluster2ce = {}
for ce in ceClusters:
    if True: #len(ceClusters[ce]) > 1:
        #print ce, ceClusters[ce]
        cluster = tuple(sorted(ceClusters[ce]))
        clusters.add(cluster)
        if cluster not in cluster2ce:
            cluster2ce[cluster] = set()
        cluster2ce[cluster].add(ce)
findSubsets = False
if findSubsets:
    keep = []
    clusters = list(clusters)
    for ii in range(len(clusters)):
        ci = set(clusters[ii])
        keepI = True
        for jj in range(len(clusters)):
            if ii != jj:
                cj = set(clusters[jj])
                keepI = keepI and not ci.issubset(cj)
        if keepI:
            keep.append(clusters[ii])
    clusters = set(keep)
for cluster in clusters:
    print len(cluster),cluster,cluster2ce[cluster]
coveredObjects =set()
for cluster in clusters:
    coveredObjects = set(cluster).union(coveredObjects)
for o in objects:
    if o not in coveredObjects:
        print o

print len(clusters)
#plt.scatter(np.array(p), np.array(npi))
#plt.show()


'''
for o1 in counts:
    for o2 in counts:
        if o1 != o2:
            for potentialCause in counts[o1]:
                if potentialCause in counts[o2]:
                    print o1,o2,potentialCause, counts[o1][potentialCause],counts[o2][potentialCause]

'''




'''
            if dat[0] == dat[4]:
                s = 'A'
            if dat[1] == dat[4]:
                s = 'B'
            if 'A' not in counts[oA]:
                counts[oA] = {}
            acause = (oB,A2B)
            aeffect = (t,s,ta,dv)
            if 'A' not in counts[oA]:
                counts[oA]['A'] = {}
            if acause not in counts[oA]:
                counts[oA]['A'][acause] = {}
            
            if aeffect not in counts[oA]['A'][acause]:
                counts[oA]['A'][acause][aeffect] = 0
            counts[oA]['A'][acause][aeffect] += 1
            
            bcause = (oA,A2B)
            beffect = (t,s,ta,dv)
            if 'B' not in counts[oB]:
                counts[oB]['B'] = {}
            if bcause not in counts[oB]['B']:
                counts[oB]['B'][bcause] = {}
            if beffect not in counts[oB]['B'][bcause]:
                counts[oB]['B'][bcause][beffect] = 0
            counts[oB]['B'][bcause][beffect] += 1
for obj1 in counts:
    for obj2 in counts:
        if obj1 != obj2:
            if 'A' in  counts[obj1] and 'A' in counts[obj2]:
                for c1 in counts[obj1]['A']:
                    for c2 in counts[obj2]['A']:
                        for e1 in counts[obj1]['A'][c1]:
                            for e2 in counts[obj2]['A'][c2]:
                                print c1,c2, counts[obj1]['A'][c1],counts[obj2]['A'][c2]


            '''
