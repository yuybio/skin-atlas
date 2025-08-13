/data/share/software/cellranger/cellranger-6.1.2/cellranger count --id=$1 \
                   --transcriptome=/data/share/ref/refdata-gex-GRCh38-2020-A \
                   --fastqs=$2 \
                   --sample=$3 \
                   --include-introns \
                   --localcores=40 \
                   --localmem=128
