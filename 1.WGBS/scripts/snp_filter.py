# coding=utf-8
#对snps.txt文件进行筛选，得到所有负链结果
from os import write
import sys
input_file=sys.argv[1]
out_file=sys.argv[2]
data=[]
seq="-1"
with open(input_file,'r')as f:
    for line in f.readlines():
        line1=line.strip('\n')
        line=line1.split('\t')
        if line[11]==seq or line[10]==seq:
            data.append(line1)
with open(out_file,'w')as o:
    o.write('\n'.join(data))
