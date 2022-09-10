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
fasterq-dump SRR11955379 --split-spot
```
# 运行cellranger
