# 从NCBI上下载初始的sra文件
以SRR11955379为例
# 将sra文件转化为fastq文件
有多种方法，fastq-dump、fasterq-dump、parallel-fastq-dump等
fastq-dump是NCBI的sra-toolkit中的工具，可将sra文件转化为fastq文件
```
fastq-dump SRR11955379 --gzip --split-files -O out
```
此处fastq-dump将结果以gzip形式输出到out文件夹中

--split-files:将双端测序分为两份,放在不同的文件,但是对于一方有而一方没有的reads直接丢弃

--split-spot:将双端测序分为两份,但是都放在同一个文件中

--split-3:将双端测序分为两份,放在不同的文件,但是对于一方有而一方没有的reads会单独放在一个文件夹里

fasterq-dump是fastq-dump的改良版，比fastq-dump速度更快，使用方式相似，但是没有gzip选项，如果想输出压缩文件，常需要与多线程压缩工具pigz联合使用

parallel-fastq-dump也可多线程运行多个fastq-dump，threads参数可以指定运行的线程数，与fastq-dump相同也可以指定压缩格式

此处使用fasterq-dump将sra文件转化为fastq文件
```
fasterq-dump SRR11955379 --split-files
```
得到SRR11955379_1.fastq和SRR11955399_2.fastq
# 将fastq文件压缩并重命名
用pigz工具设置多线程快速压缩fastq文件
```
pigz -p 12 SRR11955379_1.fastq
pigz -p 12 SRR11955379_2.fastq
```
得到SRR11955379_1.fastq.gz与SRR11955379_2.fastq.gz

为了让cellranger可以识别出得到的fastq文件，需要对fastq文件进行重命名，使其符合illumina下fastq文件的命名规则

一般命名格式为：SampleName_S1_L001_R1_001.fastq.gz

第一部分：SampleName，样本名，与上机时在Sample Sheet中填写的一致

第二部分：S1，S***，S后跟的数字与样本在Sample Sheet中的顺序一致，从1开始。不能分配到确定样本的read会归到S0（Undetermined_S0）

第三部分：L00*，泳道lane的编号

第四部分：R*，R1表示read1，R2表示read2。R1和R2为paired end reads。同一个样本的配对的FASTQ，只有这个地方不同

第五部分：001，通常为001

此处将文件重命名为Sample1_S1_L001_R1_001.fastq.gz与Sample1_S1_L001_R2_001.fastq.gz
# 运行cellranger
```
cellranger count --id=CellrangerResult --fastqs=/media/Rome/home/guest/chens/lessons/lesson2 --transcriptome=/media/Rome/home/guest/chens/refdata-gex-GRCh38-2020-A
```
id指示结果输出的文件名

fastqs指示用于分析的fastq文件的所在目录

transcriptome指示参考基因组文件

