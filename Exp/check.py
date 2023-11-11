import numpy as np
import check_strngs 
print(check_strngs.__doc__)

print(type(check_strngs.check_strngs))
a="Upuli"
b="is pretty"
d=np.array("abcdwfghijklmnpq",dtype="c")
d=np.array(b"abcdwfghijklmnpq")
print(d)
check_strngs.check_strngs(a,b,d)
print("in test.py")
print(d)




