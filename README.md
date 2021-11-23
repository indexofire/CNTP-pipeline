# CNTP-pipeline

```
# 本地CNTP菌株HZ20-124下机数据集
HZ20-124_R1.fastq.gz
HZ20-124_R2.fastq.gz

# 公共数据库数据下载
$ cat list
ERR1485221
ERR1485223
ERR1485263
ERR1485265
ERR1485314
ERR1485318
ERR3039943
ERR3039944
ERR3039945
ERR3039946
ERR3039947
ERR3039948
ERR3039949
ERR3077564
ERR3077565
ERR3077566
# 下载公共数据库数据
$ for i in `cat list`; do fastq-dump --gzip --split-3 $i; done
```

```bash
# fastqc数据指控
$ fastqc HZ20-124_R*.fastq.gz
# 查看测序结果
$ chrome 

# 基因组组装
$ shovill --R1 HZ20-124_R1.fastq.gz --R2 HZ20-124_R2.fastq.gz --trim --outdir shovill/HZ20-124
$ cp shovill/HZ20-124/contigs.fa HZ20-124.fna

# 基因组注释
$ prokka --compliant --outdir prokka/HZ20-124 --prefix HZ20-124

# insilico MLST分型
$ mlst --scheme vcholerae HZ20-124.fna

# 毒力基因扫描
$ abricate --db vfdb *.fna > vfdb.result

# 进化树构建
$ ls *.fna | parallel -k --eta snippy --ref N16961.fasta --ctgs {} --outdir snippy/{/.}
$ snippy-core --ref N16961.fasta snippy/*
$ raxmlHPC-PTREADS-AVX32 -f a -x 12345 -p 12345 -m GTRGAMMA -#1000 -s core.aln -n snps
```


```R
> library(ggtree)
> d <- read.table("genodata.txt", header=T, check.names=F)
> t <- read.tree("RAxML_bipartitions.snps")
> p1 <- ggtree(t, layout="circular") + geom_tiplab(align=T, offset=0.02) + geom_tippoint()
> 

```
