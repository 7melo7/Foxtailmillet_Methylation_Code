#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_matrix_na80_with_counter.py
仅用标准库实现：
  1. 保留“至少一组非 NA 比例 ≥80 %”的行；
  2. 每处理 10 000 行打印一次进度。
用法：
    python filter_matrix_na80_with_counter.py input.txt output.txt
"""

import sys
from collections import defaultdict

THRESHOLD = 0.6      # 80 %
NA_STRING = 'NA'     # 可根据实际改为 '.' 等
REPORT_EVERY = 500000 # 进度报告间隔

def read_header(line, sep='\t'):
    parts = line.rstrip('\n').split(sep)
    if len(parts) < 3:
        raise ValueError('列数不足3列')
    group2cols = defaultdict(list)
    for idx, col_name in enumerate(parts[2:], 2):
        group2cols[col_name].append(idx)
    return parts[0], parts[1], dict(group2cols)

def pass_filter(parts, group2cols, threshold=THRESHOLD, na=NA_STRING):
    pass_time = 0
    for grp, cols in group2cols.items():
        total = len(cols)
        non_na = sum(1 for c in cols if parts[c] != na)
        if non_na / total >= threshold:
            pass_time += 1
    return pass_time

def main():
    if len(sys.argv) != 3:
        print('Usage: python filter_matrix_na80_with_counter.py input.txt output.prefix', file=sys.stderr)
        sys.exit(1)
    in_file, out_prefix = sys.argv[1], sys.argv[2]
    
    fixed_file = out_prefix + ".fixed.txt"
    onegroup_file = out_prefix + ".onegroup.txt"

    with open(in_file, 'r', encoding='utf-8') as fin, \
         open(fixed_file, 'w', encoding='utf-8') as fout1, \
         open(onegroup_file, 'w', encoding='utf-8') as fout2:

        header = fin.readline()
        if not header:
            raise ValueError('空文件')
        fout1.write(header)
        fout2.write(header)

        chrom_col, pos_col, group2cols = read_header(header)
        kept_fixed = 0
        kept_onegroup = 0
        total = 0

        for line in fin:
            parts = line.rstrip('\n').split('\t')
            total += 1
            if total % REPORT_EVERY == 0:
                print(f'已处理 {total} 行 …')
            pass_time = pass_filter(parts, group2cols)
            if pass_time == 2:
                fout1.write(line)
                kept_fixed +=1
                continue
            if pass_time == 1:
                fout2.write(line.replace("NA", '0'))
                kept_onegroup += 1


    print(f'完成！共 {total} 行，保留 {kept_fixed} 行同时存在两个群体中 ({kept_fixed/total*100:.2f}%)')
    print(f'完成！共 {total} 行，保留 {kept_onegroup} 行只存在单个群体中 ({kept_onegroup/total*100:.2f}%)')

if __name__ == '__main__':
    main()

