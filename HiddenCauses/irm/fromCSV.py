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

collisions = set()
effects = set()
a2col = {}
b2col = {}
d2col = {}
c2eff = {}
e2s = {}
e2t = {}
e2ta = {}
e2dv = {}

def addTo(dic,k,v,val=1):
    if k not in dic:
        dic[k] = {}
    if v not in dic[k]:
        dic[k][v] = 0
    dic[k][v] += val
noEffect = set()
effects = set()
for dat in data:
    oA = dat[0]
    oB = dat[1]
    A2B = dat[2]
    t = dat[3]
    ind = 2#'OTHER'
    if dat[0] == dat[4]:
        ind = 0#'A'
    if dat[1] == dat[4]:
        ind = 1#'B'
    s = ind #dat[4]
    #ta = dat[5]
    dv = dat[6]
    
    if t != 'None' and oA != oB:
        coll = (oA,oB,A2B)
        collisions.add(coll)
        eff = (t,s,dv)
        #eff = (t,s,ta,dv)
        effects.add(eff)
        addTo(a2col,oA,coll)
        addTo(b2col,oB,coll)
        addTo(d2col,A2B,coll)
        addTo(c2eff,coll,eff)
        addTo(e2s,eff,s)
        addTo(e2t,eff,t)
        #addTo(e2ta,eff,ta)
        addTo(e2dv,eff,dv)
        effects.add(eff)
        noEffect.add((oA,oB,A2B))
    if t == 'None' and oA != oB:
        noEffect.add((oA,oB,A2B))
for pc in noEffect:
    for e in effects:
        addTo(c2eff,pc,eff,0)
    
c2i = {d:i for d,i in zip(collisions,range(len(collisions)))}
e2i = {d:i for d,i in zip(effects,range(len(effects)))}

print sys.argv[1]
print '5 1'
print len(d2ind['ObjectA']),len(d2ind['ObjectA']),1,1
print len(d2ind['A2BDir']),len(d2ind['A2BDir']),len(d2ind['A2BDir']),0
print len(d2ind['VelChange']),len(d2ind['VelChange']),len(d2ind['VelChange']),0
print 3,3,3,0
#print len(c2i),len(c2i),len(c2i),0
print len(d2ind['EffectType']),len(d2ind['EffectType']),len(d2ind['EffectType']),0
#print len(e2i),len(e2i),len(e2i),0
print '\n' # (A,B,Dir) -> (ET,Ind,V)
print 6,0,0,1,4,3,2
print '\n'

for c in c2eff:
    for eff in c2eff[c]:
        count = c2eff[c][eff]
        #if count > 0:
        #    count = count
        print 0,d2ind['ObjectA'][c[0]],d2ind['ObjectA'][c[1]],d2ind['A2BDir'][c[2]],d2ind['EffectType'][eff[0]],eff[1],d2ind['VelChange'][eff[2]]   ,   count
      


pickle.dump(d2ind,open('{}.pkl'.format(sys.argv[1]),'w'))
