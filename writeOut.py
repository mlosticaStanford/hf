import numpy as np

distances = np.arange(0.1,10.1,0.1).tolist()

finalDist, energies = [],[]
for dist in distances:
    d = '%.1f' % dist
    f = open(d + '.out', 'r')
    lines = f.readlines()
    for line in lines:
        l = line.split()
        if len(l) > 3:
            if l[0]=="Final" or l[0]=="Total":
                finalDist.append(float(d))
                energies.append(float(l[5]))
    f.close()

f = open('allE.txt', 'w')
for i in range(len(finalDist)):
    f.write('%.1f      %f\n' % (finalDist[i], energies[i]))
f.close()
