nextflow run dovetail-genomics/epiPipe-nf \
    --design ./design.csv \
    --outDir tmp1 \
    --test \
    -profile vip \
    -latest \
    -r main