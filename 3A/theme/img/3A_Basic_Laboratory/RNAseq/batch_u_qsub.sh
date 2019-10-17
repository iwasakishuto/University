#!/bin/sh
#PBS -q 1
#PBS -l select=1:ncpus=1:mpiprocs=1:ompthreads=8
#PBS -l walltime=0:30:00
#PBS -W group_list=gj29
DB_DIR="db"
RNA_DIR="RNAseq"
function DecompressHandler() {
  : '
  @params ${1} Extension 1
  @params ${2} Extension 2
  @params ${3} File Name (hoge.${1}.${2})
  '
  n_compressed_fn=${#3}
  if [ ${2} = "zip" ]; then
    unzip ${3}
    n_extension=4
  elif [ ${2} = "tar" ]; then
    tar ${3}
    n_extension=4
  elif [ ${2} = "gz" ]; then
    gunzip ${3}
    n_extension=3
  elif [ ${2} = "bz2" ]; then
    bzip2 -d ${3}
    n_extension=4
  elif [ ${2} = "lha" -o ${2} = "lzh"]; then
    lha x ${3}
    n_extension=4
  elif [ ${1} = "tar" ]; then
    tar ${3}
    n_extension=$((${#2}+5))
  else
    n_extension=0
  fi
  de_compressed_fn=${3:0:$(($n_compressed_fn-$n_extension))}
  mv $de_compressed_fn "../${DB_DIR}/"
  if [ ${n_extension} -ne 0 ]; then
    rm ${3}
  fi
  echo "${de_compressed_fn}"
}
# Download, Decompress, Dispose
function D3() {
  : '
  @params ${1} URL
  '
  fn=`wget -nv --content-disposition $1 2>&1 |cut -d\" -f2`
  extensions=( `echo $fn | tr -s '.' ' '`)
  n_extensions=${#extensions[@]}
  ext1=${extensions-2}
  ext2=${extensions-1}
  DecompressHandler $ext1 $ext2 $fn
}
#=== START ===
cd $PBS_O_WORKDIR/$RNA_DIR}
if [ ! -d $DB_DIR ]; then
  mkdir $DB_DIR
fi
if [ ! -d $RNA_DIR ]; then
  mkdir $RNA_DIR
fi
# 1.データの取得
SRA_FILE=`D3 ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/ERX/ERX288/ERX288491/ERR315326/ERR315326.sra`
mv "../${DB_DIR}/$SRA_FILE" .
fasterq-dump $SRA_FILE -v --threads 8 --split-files -O ./
# 2.品質チェック
fastqc -t 8 "${SRA_FILE}_1.fastq" "${SRA_FILE}_2.fastq"
# 3.マッピング
REF_GENOME_FILE=`D3 https://cloud.biohpc.swmed.edu/index.php/s/grch38/download`
time hisat2 -x "../${DB_DIR}/${REF_GENOME_FILE}/genome" -1 "${SRA_FILE}_1.fastq" -2 "${SRA_FILE}_2.fastq" -p 8 -S "hisat_output_ERR315326.sam"
# 4. IGVに必要なインデックスファイル作成
samtools view --threads 8 -b "hisat_output_ERR315326.sam" -o "hisat_output_ERR315326.bam"
samtools sort --threads 8 "hisat_output_ERR315326.bam" -o "hisat_output_ERR315326.sorted.bam"
# 6. リード数のカウント
GENE_ANNO_FILE=`D3 ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz`
featureCounts hisat_output_ERR315326.bam -p -t exon -g gene_id -s 0 -T 8 -BC -a "../${DB_DIR}/${GENE_ANNO_FILE}" -o "Counts_BC_ERR315326.txt"
featureCounts hisat_output_ERR315326.bam -p -t exon -g gene_id -s 0 -T 8 -MOBC -a "../${DB_DIR}/${GENE_ANNO_FILE}" -o "Counts_MOBC_ERR315326.txt"
