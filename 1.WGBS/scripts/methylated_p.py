# coding=utf-8
#针对bismark结果增加一列二项分布的p值

#file1:转化率的文件里面只有一个数值
#file2:bismark结果文件 
import sys
import os
from scipy.stats import binomtest

file1=sys.argv[1]
file2=sys.argv[2]
out_file=sys.argv[3]
cytosine_num=int(sys.argv[4])
file=open(file1)
rate=file.read().splitlines()
rate=float(''.join(rate))
file.close()
data=[]
with open(file2,'r')as f:
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        mc=int(line[3])
        c=int(line[4])
        total=mc+c
        if total >= cytosine_num:
            p=str(binomtest(mc,n=total,p=rate,alternative='greater').pvalue)
            data.append("%s\t%s"%('\t'.join(line),p))
        else:
            data.append("%s\t%s"%('\t'.join(line),'1'))
        
with open(out_file,'w')as o:
    o.write('\n'.join(data))

