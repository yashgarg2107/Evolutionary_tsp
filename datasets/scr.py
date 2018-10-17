import re

inp = open("temp","r")
out = open("ex15.txt","w+")

for line in inp:

	line = re.sub(' +', ' ',line)
	line = line.strip().split(' ')
	x = float(line[1]) 
	y = float(line[2])
	out.write(str(x)+' '+str(y)+'\n')