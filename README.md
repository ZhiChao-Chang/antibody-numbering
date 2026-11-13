# Antibody Numbering

本地离线的抗体可变区批量编号工具。输入FASTA，一次输出Kabat、Chothia、IMGT、Martin、AHo五种编号。算法由本项目实现，使用公开参考框架、序列比对和显式编号规则。

## 快速开始

需要Python 3.10或更高版本；使用标准库，无需安装第三方包、模型权重或GPU。

```bash
git clone https://github.com/ZhiChao-Chang/antibody-numbering.git
cd antibody-numbering

python3 -S scripts/number_antibodies_independent.py \
  --input examples/6vja_fv.fasta \
  --output results/example_001
```

编号自己的序列时，把`--input`替换为实际FASTA路径。每条记录应为一个抗体可变域，VH和VL分成独立记录；输入名称不影响编号。完整重/轻链建议先提取可变区。普通靶蛋白不适用抗体编号，多域串联不在当前支持范围内。

## 命令行参数

| 参数 | 说明 |
| --- | --- |
| `--input` | 必填；包含一条或多条序列的FASTA文件 |
| `--output` | 必填；尚未存在的输出目录，程序自动创建 |
| `--frameworks` | 可选；默认使用项目内的公开框架JSON |
| `--frequency-cache` | 可选；参考频率目录，可重复传入，按参数顺序查找 |

默认框架为`data/reference/independent_numbering_public_20261109/framework_templates.json`，相对于入口脚本所在项目定位。可以从其他工作目录使用入口脚本的绝对路径调用。自定义参考通过`--frameworks /path/to/framework_templates.json`显式覆盖，文件须满足现有数据格式。

需要氨基酸参考频率时：

```bash
python3 -S scripts/number_antibodies_independent.py \
  --input examples/6vja_fv.fasta \
  --output results/example_with_frequencies \
  --frequency-cache data/reference/abysis_frequencies_all_20261105 \
  --frequency-cache data/reference/independent_numbering_frequencies_20261109
```

频率是该标准位置某种氨基酸的公开参考出现比例，不是编号正确概率。缓存是部分位置的固定快照；没有对应统计时留空并记录状态，不影响编号，也不会自动联网补取。不指定频率目录时，频率列留空。

## 输出文件

| 文件 | 内容 |
| --- | --- |
| `numbering_with_frequencies.csv` | 逐序列、逐体系、逐残基的编号；包含原始位置、链类型和状态 |
| `records_summary.csv` | 每条记录、每种体系的处理状态及编号残基数 |
| `backend_outputs.json` | 详细映射、分段、框架选择和兼容说明 |
| `input_ledger.json` | 原始输入、清理后序列和位置对应关系 |
| `manifest.json` | 文件哈希、使用的框架及频率来源、运行统计 |
| `missing_frequency_positions.json` | 缺少参考频率的位置 |

`raw_index`与`normalized_index`从1开始。标准编号以字符串保存，例如`H35A`、`H112.1`。清理时移除`-`，将`?`表示为`X`；含`*`、非法字符或空序列的记录跳过并记录原因。`unclassified`表示无法识别；`unsupported`表示相关体系超出已实现规则；`numbered`不等于经过独立正确性验证，域外或部分兼容输出仍可能留空位置。

## 目录结构

```text
scripts/number_antibodies_independent.py  命令行入口
src/antibody_numbering/engine.py          核心定位、比对、编号规则
data/reference/                          默认框架与可选参考频率
examples/6vja_fv.fasta                    公开VH/VL示例
tests/                                  单元测试与公开参考用例
docs/algorithm.md                        算法结构与支持范围
docs/provenance.json                     导入文件的初始SHA-256记录
```

## 测试

```bash
python3 -B -S -m unittest discover -s tests -v
```

测试覆盖动态规划比对、环区插入、κ/λ差异、端部截短、未知残基、离线推理、负对照及命令行默认参数和覆盖行为。示例和测试数据来自公开6VJA结构；测试不是对所有抗体正确率的保证。

## 方法与限制

推理不调用ANARCI、ANARCII、HMMER或逐序列网站接口，无需训练。默认框架包含五条公开重/轻链参考，CDR序列和对应标签已遮蔽。极端环长、罕见框架插入缺失、非常规物种或严重变异可能无法识别或误判。部分旧版abYsis兼容行为被明确保留，详见[算法说明](docs/algorithm.md)。

本仓库只包含编号工具及公开参考资料；运行生成的输入副本和结果默认不纳入Git。公开资料的来源及适用边界见[数据说明](data/reference/README.md)。

## 许可

自编源代码采用[MIT许可证](LICENSE)。公开结构、编号标签及频率资料保留其来源与原权利，不因代码许可证而重新授权，详见[参考资料说明](data/reference/README.md)。
