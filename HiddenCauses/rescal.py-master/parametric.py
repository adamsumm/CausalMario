


from subprocess import call
import os
for rank in [5,8,10,20,30,50,70,100]:
    for thresh in [0,5,10,20,50,100]:
        print "python relationalize.py Mario/t_.csv {} {} 5".format(rank,thresh)
        os.system("python relationalize.py Mario/t_.csv {} {} 5".format(rank,thresh))
