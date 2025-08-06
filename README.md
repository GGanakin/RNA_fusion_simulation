# Fusion Sequence Simulator (sim_fusion_v1.py)

一个用于生成融合基因相关测序数据的Python脚本，可以模拟span、normal和split三种类型的paired-end reads。

## 功能特性

- 生成跨越融合断点的span reads
- 生成正常区域的normal reads  
- 生成跨断点分割的split reads
- 支持独立控制三种reads类型的数量
- 模拟真实的Illumina测序数据格式
- 自动引入测序错误和质量分数
- 输出压缩的FASTQ文件

## 使用方法

### 基本语法
```bash
python sim_fusion_v1.py -f <FASTA文件> -n <输出前缀> -b <断点位置> [选项]
```

### 必需参数
- `-f, --fasta`: 输入的FASTA文件路径
- `-n, --name`: 输出文件前缀
- `-b, --breakpoint`: 融合断点位置

### 可选参数
- `--span-reads`: span reads数量 (默认: 2000)
- `--normal-reads`: normal reads数量 (默认: 2000)
- `--split-reads`: split reads数量 (默认: 2000)

### 使用示例

```bash
# 使用默认reads数量
python sim_fusion_v1.py -f fusion_seq.fasta -n sample1 -b 5000

# 自定义各类reads数量
python sim_fusion_v1.py -f fusion_seq.fasta -n sample2 -b 5000 \
    --span-reads 1500 --normal-reads 3000 --split-reads 800
```

## 输出文件

脚本会生成以下文件：

### 中间文件
- `{prefix}_span_R1.fastq` / `{prefix}_span_R2.fastq` - span reads
- `{prefix}_normal_R1.fastq` / `{prefix}_normal_R2.fastq` - normal reads
- `{prefix}_split_R1.fastq` / `{prefix}_split_R2.fastq` - split reads

### 最终输出
- `{prefix}_fusion_R1.fastq.gz` - 合并压缩的R1文件
- `{prefix}_fusion_R2.fastq.gz` - 合并压缩的R2文件

## Reads类型说明

### Span Reads
- 跨越融合断点的reads对
- 模拟融合基因特征性的spanning reads
- 用于检测融合事件

### Normal Reads  
- 来自正常基因区域的reads对
- 不跨越断点，模拟背景噪音
- 按断点位置比例分配到融合基因的两侧

### Split Reads
- 部分序列跨越断点的reads
- 模拟split-read检测算法的输入
- 用于精确定位断点位置

## 技术参数

- **Read长度**: 150bp
- **插入片段大小**: 100-370bp (随机)
- **错误率**: 0.1% (0.001)
- **质量分数**: Q30-Q40
- **Read命名**: 模拟Illumina格式

## 依赖要求

- Python 3.x
- 标准库: `random`, `subprocess`, `argparse`
- 系统工具: `cat`, `gzip`

## 注意事项

1. 输入FASTA文件应包含完整的融合序列
2. 断点位置应在序列范围内
3. 确保有足够磁盘空间存储输出文件
4. 脚本会自动清理中间文件，只保留最终压缩文件

## 示例工作流程

```bash
# 1. 准备融合序列文件
# fusion_seq.fasta 包含完整的融合基因序列

# 2. 运行模拟
python sim_fusion_v1.py -f fusion_seq.fasta -n test_sample -b 8000 \
    --span-reads 2000 --normal-reads 5000 --split-reads 1000

# 3. 检查输出
ls -la test_sample_fusion_*.fastq.gz
```

这样生成的测序数据可用于测试融合基因检测算法的性能和准确性。