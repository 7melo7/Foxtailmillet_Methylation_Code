# coding=utf-8
#针对bismark结果有p和q值根据标准增加三列（甲基化mC，未甲基化C，待删除N）

import sys
file1=sys.argv[1]
file2=sys.argv[2]
data=[]
with open(file1,'r')as f:
    for line in f.readlines():
        line=line.strip('\n').split('\t')
        mc=int(line[3])
        c=int(line[4])
        total=mc+c
        q=float(line[8])
        if total<4:
            a="NA"
        elif total>=4 and q<=0.05:
            a="mC"
        elif total>=4 and q>0.05:
            a="C"
        data.append("%s\t%s"%('\t'.join(line),a))
with open(file2,'w')as o:
    o.write('\n'.join(data))
