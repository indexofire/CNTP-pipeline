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
> df <- read.table("genodata.txt", header=T, check.names=F)
> t <- read.tree("RAxML_bipartitions.snps")
> gp <- split(row.names(df), df$clade)
# 整合clade信息到树中
> t <- groupOTU(t, gp)
> p1 <- ggtree(t, layout="fan", open.angle=86, aes(color=group)) + geom_tiplab(align=T, offset=0.02) + geom_tippoint() + geom_treescale(x=0,y=0,width=0.05)
> p2 <- gheatmap(p1, df$type, width=0.1, offset=0.2, colnames_offset_y=0) + scale_fill_viridis_d(name="type")
> library(ggnewscale)
> p2 <- p2 + new_scale_fill()
> p2 <- gheatmap(p2, df[8:9], width=0.2, offset=0.3, colnames_offset_y=0) + scale_fill_viridis_d(name="gene")
> rotate_tree(p1, 88)
```
