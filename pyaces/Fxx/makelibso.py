import os 
list_dir = [name for name in os.listdir(".") if os.path.isdir(name)]

for d in list_dir:
    if (d=="Makefiles"): continue
    if (d[0]=="."): continue
    os.system("mv "+ d+"/*.a d+/*.so")
    #print("rm "+ d+"/GNU*")
