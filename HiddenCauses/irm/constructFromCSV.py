import sys
import pickle

with open(sys.argv[1],'r') as openfile:
    types = openfile.readline().rstrip().split(',')
    ind2type = {i:t for i,t in zip(range(len(types)),types)}
    typeSets = {t:set() for t in types}
    data = []
    for line in openfile:
        dat = line.rstrip().split(',')
        data.append(dat)
        for ii in range(len(dat)):
            typeSets[ind2type[ii]].add(dat[ii])
            
linked = [['ObjectA','ObjectB','Source','Target']]

for links in linked:
    urSet = set()
    for link in links:
        urSet = set(list(urSet) + list(typeSets[link]))
    for link in links:
        typeSets[link] = urSet

d2ind = {t:{d:i for d,i in zip(typeSets[t],range(len(typeSets[t])))} for t in typeSets}


# O,O,D -T-> O,O,dv 

counts = {}

def add2Counts(rel,num):
    if rel not in counts:
        counts[rel] = num
    else:
        counts[rel] += num
noEffects = set()
effects= set()
for dat in data:
    oA = dat[0]
    oB = dat[1]
    A2B = dat[2]
    t = dat[3]
    s = dat[4]
    ta = dat[5]
    dv = dat[6]

    rel = (t,oA,oB,A2B,s,ta,dv)
    #If we observe an interaction and nothing happens, we want to set all possible effects to 0

    if 'None' in t:
        noEffects.add((oA,oB,A2B))
        pass
    else:
        noEffects.add((oA,oB,A2B))
        effects.add((t,s,ta,dv))
        add2Counts(rel,1)
        pass

for c in noEffects:
    for e in effects:
        if d2ind['ObjectA'][c[0]] != d2ind['ObjectA'][e[2]] and  d2ind['ObjectA'][c[1]] != d2ind['ObjectA'][e[2]] and d2ind['ObjectA'][c[0]] != d2ind['ObjectA'][e[1]] and  d2ind['ObjectA'][c[1]] != d2ind['ObjectA'][e[1]]:
            rel = (e[0],c[0],c[1],c[2],e[1],e[2],e[3])
            add2Counts(rel,0)
print sys.argv[1]
print '4 5'
print len(d2ind['ObjectA']),len(d2ind['ObjectA']),1,1
print len(d2ind['A2BDir']),len(d2ind['A2BDir']),len(d2ind['A2BDir']),0
print len(d2ind['VelChange']),len(d2ind['VelChange']),len(d2ind['VelChange']),0
print 3,3,3,0 #A, B, other

#VelChange
print 5,0,0,1,3,2
#Add 
print 4,0,0,1,0
#Change
print 5,0,0,1,3,0
#Delete 
print 4,0,0,1,3
totalCounts = 0
for relation in counts:
    count = counts[relation]
    totalCounts += count
meanCount = float(totalCounts)/float(len(counts))
for relation in counts:
    relationType = relation[0]
    count = counts[relation]
    #print relation,count

    if  d2ind['ObjectA'][relation[1]] !=  d2ind['ObjectA'][relation[2]]:
        relation = list(relation)[1:]
        a = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[0]]
        b = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[1]]
        ind = 0
        if a:
            ind = 1
        if b:
            ind = 2
            
        if relationType == 'VelChange':
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], ind,d2ind['VelChange'][relation[5]], count
        if relationType == 'Add':
            if d2ind['ObjectA'][relation[0]] !=  d2ind['ObjectA'][relation[3]] and  d2ind['ObjectA'][relation[1]]  !=  d2ind['ObjectA'][relation[3]]:
                print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], d2ind['ObjectA'][relation[3]],  count
        if relationType == 'Delete':
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], ind,  count
        if relationType == 'Change':
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], ind, d2ind['ObjectA'][relation[4]],   count

'''   
print sys.argv[1]
print '4 5'
print len(d2ind['ObjectA']),len(d2ind['ObjectA']),1,1
print len(d2ind['A2BDir']),len(d2ind['A2BDir']),len(d2ind['A2BDir']),0
print len(d2ind['VelChange']),len(d2ind['VelChange']),len(d2ind['VelChange']),0
print 3,3,3,0 #A, B, other

#VelChange
print 4,0,0,1,2
#Add 
print 4,0,0,1,0
#Change From
print 4,0,0,1,3
#Delete 
print 4,0,0,1,3
#Change To
print 4,0,0,1,0

print '\n\n'

for relation in counts:
    relationType = relation[0]
    count = counts[relation]
    if  d2ind['ObjectA'][relation[1]] !=  d2ind['ObjectA'][relation[2]]:
        relation = list(relation)[1:]
        if relationType == 'VelChange':
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], d2ind['VelChange'][relation[5]], count
        if relationType == 'Add':
            if d2ind['ObjectA'][relation[0]] !=  d2ind['ObjectA'][relation[3]] and  d2ind['ObjectA'][relation[1]]  !=  d2ind['ObjectA'][relation[3]]:
                print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], d2ind['ObjectA'][relation[3]],  count
        if relationType == 'Delete':
            a = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[0]]
            b = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[1]]
            ind = 0
            if a:
                ind = 1
            if b:
                ind = 2
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], ind,  count
        if relationType == 'Change':
            a = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[0]]
            b = d2ind['ObjectA'][relation[3]] == d2ind['ObjectA'][relation[1]]
            ind = 0
            if a:
                ind = 1
            if b:
                ind = 2
            print d2ind['EffectType'][relationType]-1,    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], ind,   count
            print 4, d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]],   d2ind['ObjectA'][relation[4]], count

'''
pickle.dump(d2ind,open('{}.pkl'.format(sys.argv[1]),'w'))
