nextflow run dovetail-genomics/epiPipe-nf \
    --design ./design.csv \
    --outDir tmp2 \
    --genome s3://dovetail-public/publicData/nextflow-assets/genomes/bwa-mem2-index/hg38/hg38.fa \
    --test \
    -profile vip \
    -latest \
    -r main
