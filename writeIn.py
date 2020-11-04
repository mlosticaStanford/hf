import numpy as np

for i in range(100):
    a = 0.1 + (0.1*i)
    name = '%.1f' % a 
    f = open(name + '.in', 'w')
    f.write('$molecule\n')
    f.write('\n')
    f.write('3\n\n')
    f.write('H %f %f %f\n' % (0., ((3**.5)/3)*a, 0.))
    f.write('H %f %f %f\n' % (.5*a, -((3**.5)/6)*a, 0.))
    f.write('H %f %f %f\n' % (-.5*a, -((3**.5)/6)*a, 0.))
    f.write('\n')
    f.write('$calcParams \n')
    f.write('\n')
    f.write('basis = sto3g\n')
    f.write('type = uhf\n')
    f.write('tol = 0.00001 \n')
    f.write('alphaElecs = 2 \n')
    f.write('betaElecs = 2\n')
    f.write('maxIters = 1000\n')
    f.write('mix = 0')
    f.close()
