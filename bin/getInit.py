import sys,os

for i in os.listdir(sys.argv[1]):
    if '.tdd' in i:
        x=i.split(".tdd")[0].split("_")
        if len(x)<3:
#            print(x[0])
            print("uv_".join(x[0:-1]));
        else:
            print("uv_".join(x[0:-1]));
        exit()
