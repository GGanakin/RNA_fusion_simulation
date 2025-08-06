import random
import subprocess
import argparse

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence


parser = argparse.ArgumentParser(description="生成fusion相关span序列")
parser.add_argument('-f', '--fasta', required=True, help='输入的Fasta文件路径')
parser.add_argument('-n', '--name', required=True, help='输出文件前缀')
# parser.add_argument('-l', '--length', type=int, default=150, help='reads长度')
parser.add_argument('-b', '--breakpoint', type=int, required=True, help='断点位置')
parser.add_argument('--span-reads', type=int, default=2000, help='span reads数量')
parser.add_argument('--normal-reads', type=int, default=2000, help='normal reads数量')
parser.add_argument('--split-reads', type=int, default=2000, help='split reads数量')
# parser.add_argument('--fa', action='store_true',help="输入的序列是fasta格式")
args = parser.parse_args()

# 将整合后的序列赋值给 fusion_seq
fusion_seq = read_fasta(args.fasta)
all_len =  len(fusion_seq)
# read_length = args.length
read_length = 150
num_span_reads = args.span_reads  # span 配对读长数量
num_normal_reads = args.normal_reads  # normal 配对读长数量
num_split_reads = args.split_reads  # split 配对读长数量
error_rate = 0.001  # 0.1% 错误率
fusion_breakpoint = args.breakpoint  # 断点位置（chr10 和 chr8 拼接处）
output_span_fastq1 = args.name + "_span_R1.fastq" 
output_span_fastq2 = args.name + "_span_R2.fastq" 


output_nor_fastq1 = args.name + "_normal_R1.fastq" 
output_nor_fastq2 = args.name + "_normal_R2.fastq" 
nor_rate = fusion_breakpoint / all_len

output_split_fastq1 = args.name + "_split_R1.fastq" 
output_split_fastq2 = args.name + "_split_R2.fastq" 
def introduce_errors(seq, error_rate):
    """引入随机测序错误（替换）"""
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) if random.random() < error_rate else base for base in seq)

def generate_quality(length, min_q=30, max_q=40):
    """生成高精度 Phred 分数（Q30-Q40）"""
    return ''.join(chr(33 + random.randint(min_q, max_q)) for _ in range(length))

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in reversed(seq)])

def generate_read_name(read_num, is_read1=True):
    """模拟 Illumina read 名称格式"""
    instrument = "A01288"
    run_id = "199"
    flowcell = "HLHH5DSXC"
    lane = "1"
    tile = random.randint(1100, 1199)
    x_coord = random.randint(1000, 9999)
    y_coord = random.randint(1000, 9999)
    read_id = "1" if is_read1 else "2"
    filter_flag = "N"
    control_number = "0"
    barcode = ''.join(random.choice('ACTG') for _ in range(8)) + '+' + ''.join(random.choice('ACTG') for _ in range(8))
    return f"@{instrument}:{run_id}:{flowcell}:{lane}:{tile}:{x_coord}:{y_coord} {read_id}:{filter_flag}:{control_number}:{barcode}"

with open(output_span_fastq1, "w") as f1, open(output_span_fastq2, "w") as f2:
    for i in range(num_span_reads):
        insert_size = random.randint(100, 370)
        steps = random.randint(10, 140)
        # start1 = len(fusion_seq) - read_length - insert_size
        start1 = fusion_breakpoint - read_length - insert_size + steps
        read1 = fusion_seq[start1:start1 + read_length]
        read2 = reverse_complement(fusion_seq[start1 + insert_size:start1 + insert_size + read_length])
        # read2 = fusion_seq[start1 + insert_size:start1 + insert_size + read_length]
        read1 = introduce_errors(read1, error_rate)
        read2 = introduce_errors(read2, error_rate)
        qual1 = generate_quality(read_length)
        qual2 = generate_quality(read_length)
        # 生成 Illumina 格式 read 名称
        read_name1 = generate_read_name(i + 1, is_read1=True)
        read_name2 = generate_read_name(i + 1, is_read1=False)
        # 写入 FASTQ 文件
        f1.write(f"{read_name1}\n{read1}\n+\n{qual1}\n")
        f2.write(f"{read_name2}\n{read2}\n+\n{qual2}\n")

# 生成配对normal FASTQ 文件
with open(output_nor_fastq1, "w") as f1, open(output_nor_fastq2, "w") as f2:
    for i in range(num_normal_reads):
        insert_size = random.randint(100, 370)
        normal_seq = fusion_seq[:fusion_breakpoint] if random.random() < nor_rate else fusion_seq[fusion_breakpoint:]
        
        start1 = random.randint(0,len(normal_seq) - insert_size - read_length)
        read1 = normal_seq[start1:start1 + read_length]
        read2 = reverse_complement(normal_seq[start1 + insert_size:start1 + insert_size + read_length])
        read1 = introduce_errors(read1, error_rate)
        read2 = introduce_errors(read2, error_rate)
        qual1 = generate_quality(read_length)
        qual2 = generate_quality(read_length)
        # 生成 Illumina 格式 read 名称
        read_name1 = generate_read_name(i + 1, is_read1=True)
        read_name2 = generate_read_name(i + 1, is_read1=False)
        # 写入 FASTQ 文件
        f1.write(f"{read_name1}\n{read1}\n+\n{qual1}\n")
        f2.write(f"{read_name2}\n{read2}\n+\n{qual2}\n")

# 生成配对split FASTQ 文件
with open(output_split_fastq1, "w") as f1, open(output_split_fastq2, "w") as f2:
    for i in range(num_split_reads):
        insert_size = random.randint(100, 370)
        break_pos = random.randint(1, 50)
        break_read = "R1" if random.random() < 0.5 else "R2"
        start1 = fusion_breakpoint - break_pos - 50 if random.random() < 0.5 else fusion_breakpoint - break_pos - 50 - insert_size - 150
        # start1 = len(fusion_seq) - read_length - insert_size
        read1 = fusion_seq[start1:start1 + read_length]
        read2 = reverse_complement(fusion_seq[start1 + insert_size:start1 + insert_size + read_length])
        read1 = introduce_errors(read1, error_rate)
        read2 = introduce_errors(read2, error_rate)
        qual1 = generate_quality(read_length)
        qual2 = generate_quality(read_length)
        # 生成 Illumina 格式 read 名称
        read_name1 = generate_read_name(i + 1, is_read1=True)
        read_name2 = generate_read_name(i + 1, is_read1=False)
        # 写入 FASTQ 文件
        f1.write(f"{read_name1}\n{read1}\n+\n{qual1}\n")
        f2.write(f"{read_name2}\n{read2}\n+\n{qual2}\n")

print(f"span FASTQ 文件已生成：{output_span_fastq1}, {output_span_fastq2}")
print(f"normal FASTQ 文件已生成：{output_nor_fastq1}, {output_nor_fastq2}")
print(f"split FASTQ 文件已生成：{output_split_fastq1}, {output_split_fastq2}")

output_fastq1 = args.name + "_fusion_R1.fastq" 
output_fastq2 = args.name + "_fusion_R2.fastq" 

fastq1_files = [output_span_fastq1, output_split_fastq1, output_nor_fastq1]  # 替换为你的 fastq1 文件列表
fastq2_files = [output_span_fastq2, output_split_fastq2, output_nor_fastq2]  # 替换为你的 fastq2 文件列表

# 使用 cat 合并 fastq1 文件
with open(output_fastq1, "w") as f_out:
    subprocess.run(["cat"] + fastq1_files, stdout=f_out)

# 使用 cat 合并 fastq2 文件
with open(output_fastq2, "w") as f_out:
    subprocess.run(["cat"] + fastq2_files, stdout=f_out)

# 使用 gzip 压缩 fastq 文件
subprocess.run(["gzip", output_fastq1])
subprocess.run(["gzip", output_fastq2])

print(f"所有 fastq1 文件已合并并压缩为：{output_fastq1}.gz")
print(f"所有 fastq2 文件已合并并压缩为：{output_fastq2}.gz")


