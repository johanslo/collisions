from collisions.neighbour import *
import matplotlib.pyplot as plt
import numpy as np

dump = NeighbourList('./tmp.dump', f_index = 5)

# print(dump.allPairlists[0])
# print(dump.allPairlists[1])

i = 0
f_pair = []
t = dump.timesteps

# for pairlist in dump.allPairlists:
#     print(i)
#     i += 1
#     print(pairlist.pairlist[5])
#     f_pair.append(pairlist.pairlist[5].pairforce)
k = 25
p = dump.collisionList[k].pair
print("Order: ", dump.collisionList[k].order)
t, f_pair = dump.extractEvent(p)

for col in dump.collisionList :
    print(i, ":", col.pair, ";", col.order)
    i += 1

plt.plot(t,np.abs(f_pair))
plt.savefig("coll.png")
plt.close()
# pair = Pair(1,2,69.42,True)
# pair2 = Pair(1,3,42.69, False)
# pl = Pairlist()
# pl.addToPairlist(pair)
# pl.addToPairlist(pair2)

# print(pl)