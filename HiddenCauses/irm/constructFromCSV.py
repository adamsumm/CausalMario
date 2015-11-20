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
        for t in typeSets['EffectType']:
            if 'None' not in t:
                for s in typeSets['ObjectA']:
                    for ta in typeSets['ObjectB']:
                        for dv in typeSets['VelChange']:
                            rel = (t,oA,oB,A2B,s,ta,dv)
                            add2Counts(rel,0)
    else:
        add2Counts(rel,1)
for relation in counts:
    relationType = relation[0]
    count = counts[relation]
    
    relation = list(relation)[1:]
    print d2ind['EffectType'][relationType],    d2ind['ObjectA'][relation[0]], d2ind['ObjectA'][relation[1]],    d2ind['A2BDir'][relation[2]], d2ind['ObjectA'][relation[3]],    d2ind['ObjectA'][relation[4]], d2ind['VelChange'][relation[5]], count

pickle.dump(d2ind,open('{}.pkl'.format(sys.argv[0]),'w'))
