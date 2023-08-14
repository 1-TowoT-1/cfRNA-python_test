import sys
import re
import os

sample=sys.argv[1] #文件输入路径+样本名称
output=sys.argv[2] #文件输出结果

sample_raw_R1=f'{sample}_R1.stat'
sample_raw_R2=f'{sample}_R2.stat'
sample_clean_R1=f'{sample}_clean_R1.stat'
sample_clean_R2=f'{sample}_clean_R2.stat'

with open(sample_raw_R1,'r') as f1:
    for line in f1:
        if re.match(r'Total bases',line):
            r1_raw_base=re.match(r"Total bases : (\d+)",line).group(1)
        if re.match(r'Total reads',line):
            r1_raw_reads=re.match(r"Total reads : (\d+)",line).group(1)
        if re.match(r'% bases >=Q20',line):
            r1_raw_q20=re.match(r"% bases >=Q20 : (\d+.\d+)",line).group(1)
        if re.match(r'% bases >=Q30',line):
            r1_raw_q30=re.match(r"% bases >=Q30 : (\d+.\d+)",line).group(1)
        
with open(sample_raw_R2,'r') as f2:
    for line in f2:
        if re.match(r'Total bases',line):
            r2_raw_base=re.match(r"Total bases : (\d+)",line).group(1)
        if re.match(r'Total reads',line):
            r2_raw_reads=re.match(r"Total reads : (\d+)",line).group(1)
        if re.match(r'% bases >=Q20',line):
            r2_raw_q20=re.match(r"% bases >=Q20 : (\d+.\d+)",line).group(1)
        if re.match(r'% bases >=Q30',line):
            r2_raw_q30=re.match(r"% bases >=Q30 : (\d+.\d+)",line).group(1)

with open(sample_clean_R1,'r') as f3:
    for line in f3:
        if re.match(r'Total bases',line):
            r1_clean_base=re.match(r"Total bases : (\d+)",line).group(1)
        if re.match(r'Total reads',line):
            r1_clean_reads=re.match(r"Total reads : (\d+)",line).group(1)
        # if re.match(r'% bases >=Q20',line):
        #     r1_clean_q20=re.match(r"% bases >=Q20 : (\d+)",line).group(1)
        # if re.match(r'% bases >=Q30',line):
        #     r1_clean_q20=re.match(r"% bases >=Q30 : (\d+)",line).group(1)        


with open(sample_clean_R2,'r') as f4:
    for line in f4:
        if re.match(r'Total bases',line):
            r2_clean_base=re.match(r"Total bases : (\d+)",line).group(1)
        if re.match(r'Total reads',line):
            r2_clean_reads=re.match(r"Total reads : (\d+)",line).group(1)
        # if re.match(r'% bases >=Q20',line):
        #     r2_clean_q20=re.match(r"% bases >=Q20 : (\d+)",line).group(1)
        # if re.match(r'% bases >=Q30',line):
        #     r2_clean_q20=re.match(r"% bases >=Q30 : (\d+)",line).group(1)


all_raw_bases=int(r1_raw_base)+int(r2_raw_base)
all_raw_reads=int(r1_raw_reads)+int(r2_raw_reads)
all_clean_bases=int(r1_clean_base)+int(r2_clean_base)
all_clean_reads=int(r1_clean_reads)+int(r1_clean_reads)

raw_Q20=str(round((float(r1_raw_q20)+float(r2_raw_q20))/2,2))+"%"
raw_Q30=str(round((float(r1_raw_q30)+float(r2_raw_q30))/2,2))+"%"
sample_name=os.path.basename(sample)

with open(output,'w') as f5:
    print("Sample","raw_reads","raw_bases","clean_reads","clean_bases","Q20","Q30",sep='\t',file=f5)
    print(sample_name,all_raw_reads,all_raw_bases,all_clean_reads,all_clean_bases,raw_Q20,raw_Q30,sep='\t',file=f5)
