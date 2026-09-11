# coding=utf-8
# 统计甲基化mC，未甲基化C，待删除N
import sys
import gzip

gz_file = sys.argv[1]
ID = sys.argv[2]

c_dict = {'mC':0, 'C':0, 'NA':0}

with gzip.open(gz_file, 'rt', encoding="utf-8") as f:
    line = f.readline().rstrip()
    while line:
        if line.startswith("NC"):
            line = f.readline().rstrip()
        else:
            break
    while line:
        fields = line.split('\t')
        c_dict[fields[9]] += 1
        line = f.readline().rstrip()

with open(sys.argv[3], 'a')as f:
    f.write('\t'.join([ID, str(c_dict["NA"]), str(c_dict["C"]+c_dict["mC"]), str(c_dict["mC"]), str(round(c_dict["mC"]/(c_dict["mC"]+c_dict["C"]), 6))]) + '\n')
