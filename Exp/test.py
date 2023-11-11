import numpy as np
from numpy import array 
import add
import check_strngs 
print(add.add.__doc__)
print(add.internal.__doc__)
print(add.check_strngs.__doc__)

a=np.array([5.0,4.0,3.0])
b=np.array([3.0,4.0,5.0])

print(a)
print(b)

c=add.add(a,b,3)

print("@-test.py")
print(a)
print(b)
d="Upuli"
e="is pretty"
#g=np.array(["ijklmnpqrstupqia","ijklmnpqrstupqia"],dtype="c")
g=np.array(["ijklmnpqrstupqia","ijklmnpqrstupqia"],dtype="c")
print(g)
print(type(g))
#print(g.shape(0))
check_strngs.check_strngs(d,e,g,2)
print("in test.py")
print(g)


