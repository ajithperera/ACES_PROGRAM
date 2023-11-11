import numpy as np

#dt=np.dtype([("name",np.unicode_,16), ("grade",np.float64,(2,))])
#a["name"]
#a["grade"]

dt=np.dtype([("name",np.unicode_,16)])
a=np.array([],dtype=dt)

print(type(a))
print(a.shape)



