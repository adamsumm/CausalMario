import sys
import logging
logging.basicConfig(level=logging.INFO)
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
import numpy as np
from rescal import rescal_als
import matplotlib.pyplot as plt
import numpy as np
from numpy import dot, array, zeros, setdiff1d
from numpy.linalg import norm
from numpy.random import shuffle
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
from sklearn.metrics import precision_recall_curve, auc
import pickle
import os.path
import operator
from sklearn.preprocessing import normalize

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

def addTo(dic,k,v):
    if k not in dic:
        dic[k] = {}
    if v not in dic[k]:
        dic[k][v] = 0
    dic[k][v] += 1
    
for dat in data:
    oA = dat[0]
    oB = dat[1]
    A2B = dat[2]
    coll = (oA,oB,A2B)
    collisions.add(coll)
    t = dat[3]
    ind = 'OTHER'
    if dat[0] == dat[4]:
        ind = 'A'
    if dat[1] == dat[4]:
        ind = 'B'
    s = dat[4]
    ta = dat[5]
    dv = dat[6]
    eff = (t,s,ta,dv)
    effects.add(eff)
    addTo(a2col,oA,coll)
    addTo(b2col,oB,coll)
    addTo(d2col,A2B,coll)
    addTo(c2eff,coll,eff)
    addTo(e2s,eff,s)
    addTo(e2t,eff,t)
    addTo(e2ta,eff,ta)
    addTo(e2dv,eff,dv)
    



relations = [a2col,b2col,d2col,c2eff,e2s,e2t,e2ta,e2dv]

entities = set()
for rel in relations:
    for e in rel:
        entities.add(e)
        for t in rel[e]:
            entities.add(t)
            
e2i = {e:i for e,i in zip(entities,range(len(entities)))}


sparse = [lil_matrix((len(entities),len(entities))) for r in relations]
ind = 0
for rel in relations:
    for s in rel:
        for t in rel[s]:
            sparse[ind][e2i[s],e2i[t]] = rel[s][t]
    ind += 1 
#plt.show()

fname = '{}{}.pkl'.format(sys.argv[1],sys.argv[2])

if os.path.isfile(fname) :
    A,R,fit,itr,exectimes = pickle.load(open(fname,'r'))
else:
    A, R, fit, itr, exectimes = rescal_als(sparse, int(sys.argv[2]), init='nvecs', conv=1e-6, lambda_A=1, lambda_R=1)
    pickle.dump([A, R, fit, itr, exectimes],open(fname,'w'))

n = A.shape[0]
P = zeros((n, n, len(R)))
for k in range(len(R)):
    P[:, :, k] = dot(A, dot(R[k], A.T))
print A.shape
print len(R),R[0].shape
plt.matshow(P[:,:,0])
#plt.show()
A_ = np.mean(A,axis=0)
e2vec = {}
for e in a2col:
    e2vec[e] = A[e2i[e],:]
print e2vec
# Type 1
for e in e2vec:
    cosSimilarity = {}
    e_ = normalize(e2vec[e].reshape(-1, 1)).ravel()
    for e2 in e2vec:
        if e != e2:
            e2_ = normalize(e2vec[e2].reshape(-1, 1)).ravel()
            cosSimilarity[e2] = np.dot(e_,e2_)
        
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    #print e,sorted_x[-1],sorted_x[-2],sorted_x[-2]
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    for ii in range(len(sorted_x)-1,-1,-1):
        if sorted_x[ii][1] > float(sys.argv[4]):
            print e,sorted_x[ii]
# Type 2
membership = {}

for e in e2vec:
    #print e, e2vec[e],np.percentile(A,float(sys.argv[4]),axis=0)
    groups = (e2vec[e]>(np.percentile(A,float(sys.argv[4]),axis=0))).nonzero()#groups = (e2vec[e]>(A_*(1-float(sys.argv[4]))+A_max*float(sys.argv[4]))).nonzero()#(e2vec[e]>float(sys.argv[4])).nonzero()
    for g in groups[0]:
        if g not in membership:
            membership[g] = set()
        membership[g].add(e)
print membership
for ii in sorted(membership.keys()):
    if ii in membership:
        print ii, len(membership[ii]), membership[ii]
   
exit()     
# Type 1
groups = set()
for e in e2vec:
    cosSimilarity = {}
    e_ = normalize(e2vec[e].reshape(1, -1)).ravel()
    for e2 in e2vec:
        if e != e2:
            e2_ = normalize(e2vec[e2].reshape(1, -1)).ravel()
            cosSimilarity[e2] = np.dot(e_,e2_)
        
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    #print e,sorted_x[-1],sorted_x[-2],sorted_x[-2]
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    group = set()
    for ii in range(len(sorted_x)-1,-1,-1):
        if sorted_x[ii][1] > float(sys.argv[3]):
            group.add(sorted_x[ii][0])
    groups.add(tuple(sorted(list(group))))


keep = []
for group in groups:
    shouldKeep = True
    for oGroup in groups:
        if group != oGroup and set(group).issubset(set(oGroup)):
            shouldKeep = False
            break
    if shouldKeep:
        keep.append(group)
groups = keep
for group in groups:
    print len(group),group

            #print e,sorted_x[ii]
fname = '{}{}_{}norm.pkl'.format(sys.argv[1],sys.argv[2],sys.argv[3])

if os.path.isfile(fname) :
    A,R,fit,itr,exectimes = pickle.load(open(fname,'r'))
else:
    for ii in range(len(sparse)):
        sparse[ii] = sparse[ii] > float(sys.argv[3])
    A, R, fit, itr, exectimes = rescal_als(sparse, int(sys.argv[2]), init='nvecs', conv=1e-6, lambda_A=1, lambda_R=1)
    pickle.dump([A, R, fit, itr, exectimes],open(fname,'w'))

n = A.shape[0]
P = zeros((n, n, len(R)))
for k in range(len(R)):
    P[:, :, k] = dot(A, dot(R[k], A.T))
print A.shape
print len(R),R[0].shape
plt.matshow(P[:,:,0])
#plt.show()
A_ = np.mean(A,axis=0)
A_max = np.max(A,axis=0)
e2vec = {}
for e in a2col:
    e2vec[e] = A[e2i[e],:]

# Type 1
for e in e2vec:
    cosSimilarity = {}
    e_ = normalize(e2vec[e].reshape(-1, 1)).ravel()
    for e2 in e2vec:
        if e != e2:
            e2_ = normalize(e2vec[e2].reshape(-1, 1)).ravel()
            cosSimilarity[e2] = np.dot(e_,e2_)
        
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    #print e,sorted_x[-1],sorted_x[-2],sorted_x[-2]
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    for ii in range(len(sorted_x)-1,-1,-1):
        if sorted_x[ii][1] > float(sys.argv[4]):
            print e,sorted_x[ii]
# Type 2
membership = {}

for e in e2vec:
    #print e, e2vec[e],np.percentile(A,float(sys.argv[4]),axis=0)
    groups = (e2vec[e]>(np.percentile(A,float(sys.argv[4]),axis=0))).nonzero()#groups = (e2vec[e]>(A_*(1-float(sys.argv[4]))+A_max*float(sys.argv[4]))).nonzero()#(e2vec[e]>float(sys.argv[4])).nonzero()
    for g in groups[0]:
        if g not in membership:
            membership[g] = set()
        membership[g].add(e)
print membership
for ii in sorted(membership.keys()):
    if ii in membership:
        print ii, membership[ii]
groups = set()
for e in e2vec:
    cosSimilarity = {}
    e_ = normalize(e2vec[e].reshape(1, -1)).ravel()
    for e2 in e2vec:
        if e != e2:
            e2_ = normalize(e2vec[e2].reshape(1, -1)).ravel()
            cosSimilarity[e2] = np.dot(e_,e2_)
        
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    #print e,sorted_x[-1],sorted_x[-2],sorted_x[-2]
    sorted_x = sorted(cosSimilarity.items(), key=operator.itemgetter(1))
    group = set()
    for ii in range(len(sorted_x)-1,-1,-1):
        if sorted_x[ii][1] > float(sys.argv[4]):
            group.add(e)
    groups.add(tuple(list(group)))

for group in groups:
    print len(group),group
