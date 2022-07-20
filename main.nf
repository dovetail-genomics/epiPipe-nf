#!/usr/bin/env nextflow

params.help    = false
params.test    = false

params.design  = false
params.outDir  = false

params.genome  = 'hg38'
params.mapQ    = 40

params.resolutions   = false
params.ABresolutions = false


def helpMessage() {
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

        alignPerDesign.nf  --design ~/path/to/bam/location/design.csv --outDir ~/path/to/where/to/save/files

    Mandatory arguments:
        --outDir [path]              Path to a diectory to save the bigwig coveage files (can be local or valid S3 location.

        --design [path]              Path to a design files in csv with header has first line (id,rep,R1,R2). First col: id, second col: replicate, third col: path to R1, fourth col: path to R2
    
    Alignment:
        --genome [str]               Name or path of the genome to use. Possible choice: hg38, hg19, mm10, dm3. Default: hg38.
        --mapQ [int]                 Quqlity score to filter reads. Integer between 0 and 60. Default: 40.

    Additional parameers:
        --resolutions [integers]     Arrowhead and hiccup resolutions. Comma-seperated list of resolutions in kb for loops finding. Default: [5,10]
        --ABresolutions [integers]   ABcomp resolutions. Comma-seperated list of resolutions in kb. Default: [32,64,128]
    
    Chicago parameters:
        --resolutions [integers]     Comma-seperated list of resolutions for computing genomic bins. Default: [5,10,20]
        --panel [string]             Id of the panel to use. Valid id: hs_pc_1, mm_pc_1. Default: hs_pc_1.

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

if (!(params.outDir && params.design)) {
    exit 1, "--outDir and --design is are required arguments. Use --help to get the full usage." 
} else {
    outDir = params.outDir
}

if (params.resolutions){
    resolutions = params.resolutions.split(/,/,-1)
} else {
    resolutions = [5,10]
}

if (params.ABresolutions){
    ABresolutions = params.ABresolutions.toString().split(/,/,-1)
} else {
    ABresolutions = [32,64,128]
}

////////////////////////////////////////////////////////////////////////////////
//
// Global Variables
//
////////////////////////////////////////////////////////////////////////////////
genomeAssets = 's3://tower-ops-2/genomes/bwa-mem2/'

////////////////////////////////////////////////////////////////////////////////
//
// Log File Creation
//
////////////////////////////////////////////////////////////////////////////////
def summary = [:]
summary['Run Name']       = workflow.runName
summary['Working dir']    = workflow.workDir
summary['Genome dir']     = genomeAssets
summary['Current home']   = "$HOME"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


////////////////////////////////////////////////////////////////////////////////
//
// START OF THE PIPELINE
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//
// 1) Setup basic input channels
//
////////////////////////////////////////////////////////////////////////////////

// Creating fastq channel from design files
Channel
    .fromPath(params.design, checkIfExists: true)
	.ifEmpty { exit 1, "Can't locate desgin file: ${params.design}" }
    .splitCsv(header: ['id', 'rep', 'R1', 'R2'], skip: 1)
    .map{it -> tuple(it.id, it.rep, file(it.R1), file(it.R2))}
    .set{ fqs_ch }

if (params.genome =~ /hg19|hg38|mm10|rn6|susScr11|dm3/){
    Channel
	    .fromFilePairs("${genomeAssets}/${params.genome}/*.{0123,amb,ann,bwt.2bit.64,pac,fa}", size: -1)
        .set { bwa_index_ch }
    
    Channel
	    .fromPath("${genomeAssets}/${params.genome}/${params.genome}.fa")
	    .set { abcomp_genome_ch }

} else {
    Channel
	    .fromPath(params.genome, checkIfExists: true)
	    .ifEmpty { exit 1, "Genome not found: ${params.genome}" }
	    .set { abcomp_genome_ch }

    process  bwa_index {
        label 'major'
        container 'dovetailg/bwa-mem2'

        input:
        path(ref) from filtCont_idx_ch
        
        output:
        tuple val(id), path("*") into bwa_index_ch

        script:
        id='genomeIDX'
        """
        bwa-mem2 index -p ${id} ${ref}
        """
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// 2) Pipeline processes
//
////////////////////////////////////////////////////////////////////////////////
process bwa_mem2 {
    tag "_${prefix}"
    label 'cpu'
    container 'dovetailg/bwa-mem2'
    
    input:
    tuple val(id), val(rep), path(R1), path(R2) from fqs_ch
    tuple val(index), path(index_files) from bwa_index_ch.first() 
    
    output:
    tuple val(id), val(rep), val(prefix), path("*.bam") into  bam_chrSize_ch
    
    script:
    prefix = R1.name.toString().replaceFirst('\\..+',"")
    if (!params.test)
        """
        bwa-mem2 mem -5SP -t ${task.cpus} \
            ${index} \
            ${R1} \
            ${R2} \
        |samtools view -@ ${task.cpus} -Shb -o ${prefix}.bam - 
        """
    else 
        """
        bwa-mem2 mem -5SP -t ${task.cpus} \
    	    ${index} \
    	    <(cat ${R1}|head -n 400000) \
    	    <(cat ${R2}|head -n 400000) \
	    |samtools view -@ ${task.cpus} -Shb -o ${prefix}.bam - 
        """
}

process chr_size1 {
    tag "_${prefix}"
    label 'cpu'
    container 'dovetailg/bwa-mem2'
    
    input:
    tuple val(id), val(rep), val(prefix), path(bam) from  bam_chrSize_ch
    
    output:
    tuple val(id), val(rep), val(prefix), path(bam), path("*.tsv") into pairtools_parse_ch
    
    script:
    """
	samtools view -H ${bam} \
	|awk -v OFS='\t' '/^@SQ/ && !(\$2 ~ /:(chr|"")M/) {split(\$2,chr,":");split(\$3,ln,":");print chr[2],ln[2]}' \
	|sort -V -k1,1 > chr_size.tsv
    """
}

process pairtools_parse {
    tag "_${id}"
    label 'bastard'
    container 'dovetailg/pairtools'
    
    input:
    tuple val(id), val(rep), val(prefix), path(sam), path(chr_sizes) from pairtools_parse_ch

    output:
    tuple val(id), val(rep), val(prefix), path("*.pairsam.gz") into pairsam_part_ch

    script:
    """
    pairtools parse \
	--min-mapq ${params.mapQ} \
	--walks-policy 5unique \
	--max-inter-align-gap 30 \
	--nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--chroms-path ${chr_sizes} \
	--output ${prefix}.pairsam.gz \
	${sam} 
    """
}

process pairtools_merge_lane {
    tag "_${prefix_out}"
    label 'bastard'    
    container 'dovetailg/pairtools'
    
    input:
    tuple val(id), val(rep), val(prefix), path(sam) from pairsam_part_ch
	.groupTuple(by: [1,0])
    
    output:
    tuple val(id), val(prefix_out), path("${prefix_out}*pairsam.gz") into pairsam_ch

    script:
    prefix_out="${id}-${rep}"
    if (sam.sort().size() >1) {
	"""
	pairtools merge -o ${prefix_out}.pairsam.gz --nproc ${task.cpus} ${sam}
	"""
    } else {
	"""
	ln -sf ${sam} ${prefix_out}_ML.pairsam.gz
	"""
    }   
}

process pairtools_sort {
    tag "_${prefix}"
    label 'macro'
    container 'dovetailg/pairtools'

    input:
    tuple val(id), val(prefix), path(sam) from pairsam_ch

    output:
    tuple val(id), val(prefix), path("*_sorted.pairsam.gz") into sorted_ps_ch
    
    script:
    """
    mkdir -p tmp 
    pairtools sort --tmpdir ./tmp  \
	--nproc ${task.cpus} \
	--output ${prefix}_sorted.pairsam.gz \
	$sam 
    """
}

process pairtools_dedup {
    tag "_${prefix}"
    label 'bastard'
    container 'dovetailg/pairtools'
    
    publishDir "${params.outDir}/pairtools_stat",
    	mode: 'copy',
    	saveAs: {filename -> filename.endsWith('.stats') ? filename : null}
    
    input:
    tuple val(id), val(prefix), path(sam) from sorted_ps_ch

    output:
    tuple val(id), val(prefix), path("*_dedup.pairsam.gz") into dedup_ps_ch
    tuple val(id), val(prefix), path("*_unmapped.pairsam.gz") into unmapped_ps_ch
    path("*_pairtools.stats") into ps_stats_ch

    script:
    """
    pairtools dedup --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--mark-dups \
	--output-stats ${prefix}_pairtools.stats  \
	--output ${prefix}_dedup.pairsam.gz \
	--output-unmapped ${prefix}_unmapped.pairsam.gz \
	${sam}
    """
}


process pairtools_stats_merge {
    tag "_${id}"
    label 'minor'
    container 'dovetailg/pt-stats'

    publishDir "${params.outDir}",
	mode: 'copy'
    
    input:
    path(stats) from ps_stats_ch
	.collect()
    
    output:
    path('pairtoolsStats.csv') into merged_stats_ch
    
    script:
    """
    pairtoolsStat.sh ${stats} > pairtoolsStats.csv
    """
}

process pairtools_split_dedup {
    tag "_${prefix}"
    label 'bastard'
    container 'dovetailg/pairtools'

    input:
    tuple val(id), val(prefix), path(sam) from dedup_ps_ch
    
    output:
    tuple val(id), val(prefix), path("*.bam") into bam_parts_ch
    tuple val(id), val(prefix), path("*.valid.pairs.gz") into pairs_parts_ch, pairs_parts_ch_test

    script:
    """
    pairtools split --nproc-in ${task.cpus} --nproc-out ${task.cpus} \
	--output-sam ${prefix}_PT.bam  \
	--output-pairs ${prefix}_PT.valid.pairs.gz  \
	${sam}
    """
}

process merge_bam {
    tag "_${id}"
    label 'macro'
    container 'dovetailg/bwa-samtools'
    
    input:
    tuple val(id), val(rep), path(bam_part) from bam_parts_ch
	.groupTuple()

    output:
    tuple val(id), path("*.bam") into merged_bam_sort_ch

    script:
    bam_files = bam_part.sort()
    if (bam_files.size() >1) {
	"""
	samtools merge -@ ${task.cpus} ${id}_MB.bam ${bam_part}
	"""
    } else {
	"""
	ln -s ${bam_part} ${id}_MB.bam
	"""
    }
}

process merge_pairs {
    tag "_${id}"
    cpus 14
    memory '40 GB'
    container 'dovetailg/pairtools'
    
    publishDir "${params.outDir}/validPairs",
    	mode: 'copy'
    
    input:
    tuple val(id), val(rep), path(pairs) from pairs_parts_ch
	.groupTuple()
    
    output:
    tuple val(id), path("*.valid.pairs.gz"), path("*.px2") into pairs_chrSize_ch

    script:
    pair_files = pairs.sort()
    if (pair_files.size() >1) {
	"""
	pairtools merge -o ${id}.valid.pairs.gz --nproc ${task.cpus}  ${pairs}
	pairix ${id}.valid.pairs.gz
	"""
    } else {
	"""
	ln -s ${pairs} ${id}.valid.pairs.gz
	pairix ${id}.valid.pairs.gz
	"""
    }
}

process bam_sort {
    tag "bam_sort_${id}"
    label 'large'
    container 'dovetailg/bwa-samtools'
    
    publishDir "${params.outDir}/bam",
	mode: 'copy',
	pattern: "${id}.bam"
        
    input:
    tuple val(id), path(bam) from merged_bam_sort_ch
    
    output:
    tuple val(id), path("${id}.bam"),path("${id}.bam.bai") into bam_bigwig_ch,
	bam_capStats_ch,
	bam_cleanBam_ch,
	bam_mapFile_ch

    script:
    """
    samtools sort -m 2G \
	-@ ${task.cpus} \
	-o ${id}.bam \
	${bam} 

    samtools index -@${task.cpus} ${id}.bam
    """
}

process chr_size {
    tag "_${id}"
    label 'alto'
    container 'dovetailg/pairtools'
    
    input:
    tuple val(id), path(pairs), path(idx) from pairs_chrSize_ch
    
    output:
    tuple val(id), path(pairs), path(idx), path("*.tsv") into pairs_ch_cooler, pairs_ch_juicer
    
    script:
    """
    pairix -H -f ${pairs} \
	| awk -v OFS='\t' '/^#chromsize/  {print \$2,\$3}' \
	| sort -V -k1,1 \
	> chr_size.tsv
    """
}


process cooler_cload {
    tag "_${id}"
    label 'large'
    container 'dovetailg/cooler'

    input:
    tuple val(id), path(pairs), path(idx), path(chr_sizes) from pairs_ch_cooler
    
    output:
    tuple val(id), path("*.cool") into balance_cooler_ch
        
    script:
    """
    cooler cload pairix \
	-p ${task.cpus} \
	${chr_sizes}:1000 \
	${pairs} \
	${id}.cool
    """
}

process balance_cooler {
    tag "_${id}"
    label 'large'
    container 'dovetailg/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple val(id), path(cooler) from balance_cooler_ch

    output:
    tuple val(id), path(cooler) into zoomify_cooler_ch 
    
    script:
    """
    cooler balance --force -p ${task.cpus} ${cooler}
    """
}

process cooler_zoomify {
    tag "_${id}"
    label 'large'
    container 'dovetailg/cooler'
    
    publishDir "${params.outDir}/coolerFiles",
    	mode: 'copy'
    
    input:
    tuple val(id), path(cooler) from zoomify_cooler_ch

    output:
    tuple val(id), path("*.mcool") into mustache_mcool_ch, abcomp_mcool_ch
    
    script:
    """
    cooler zoomify --balance -p ${task.cpus} ${cooler}
    """
}

process bam2bw {
    tag "_${id}"
    label 'malform'
    container 'dovetailg/r-cov'
    
    publishDir "${params.outDir}/bigwigs",
    	mode: 'copy'
    
    input:
    tuple val(id), path(bam),path(idx) from bam_bigwig_ch
        
    output:
    tuple val(id), path ("*.bw") into bigwig_out_ch

    script:
    """
    bam2bw ${bam} ${id}.bw ${task.cpus}
    """
}

process juicer {
    tag "_${id}"
    label 'malform'
    container 'dovetailg/juicer'
    
    publishDir "${params.outDir}/hicFiles",
    mode: 'copy'
    
    input:
    tuple val(id), path(pairs), path(idx), path(chr_sizes) from pairs_ch_juicer
    
    output:
    tuple val(id), path("*.hic") into arrowhead_ch, hiccups_ch

    script:
    """
    java -Xmx96000m -Djava.awt.headless=true \
	-jar /juicer_tools.jar pre \
	--threads ${task.cpus} \
	-j ${task.cpus} \
	-k VC,VC_SQRT,KR,SCALE \
	${pairs} \
	${id}.hic \
	${chr_sizes}
    """
}

process arrowhead {
    tag "_${id}"
    label 'bastard'
    container "dovetailg/juicer"
    
	publishDir "${params.outDir}/arrowHead",
	mode: 'copy'
    
    input:
    tuple val(id), path(hic), val(res) from arrowhead_ch
	.combine(Channel.from(resolutions))
    
    output:
    tuple val(id), path("${id}_${res}kb") into arrowhead_out_ch
    
    script:
    bpRes = res.toInteger() * 1000
    """
    mkdir -p ${id}_${res}kb && touch ${id}_${res}kb/${bpRes}_blocks.bedpe
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	arrowhead \
	--threads ${task.cpus} \
	--ignore-sparsity \
	-r ${bpRes} \
	-k KR \
	${hic} \
	${id}_${res}kb
    """
}

process hiccups {
    tag "_${id}"
    label 'gpu'
    container "dovetailg/hiccups-gpu"
    
    publishDir "${params.outDir}/hiccups/",
	mode: 'copy'
    
    input:
    tuple val(id), path(hic), val(res)  from hiccups_ch
        .combine(Channel.from(resolutions.collect{it*1000}.join(',')))
    
    output:
    tuple val(id), path("${id}_loops") 
    
    script:
    """
    java -Xmx24000m \
	-jar /juicer_tools.jar \
	hiccups \
	--threads ${task.cpus} \
	    --ignore-sparsity \
	    -m 500 \
	-r ${res} \
	-k KR \
	${hic} \
	${id}_loops
    """
}

process mustache {
    tag "_${id}"
    label 'malform'
    container "dovetailg/mustache"
    
    publishDir "${params.outDir}/mustache",
	mode: 'copy'
    
    input:
    tuple val(id), path(mcool), val(res)  from mustache_mcool_ch
	.combine(Channel.from(1000,4000,16000))
    
    output:
    tuple val(id), path("*.tsv") 
    
    script:
    """
    touch ${id}_${res}kb_loops.tsv 
    mustache -p ${task.cpus} \
	-f ${mcool} \
	-r ${res} \
	-o ${id}_${res}kb_loops.tsv
    """
}

process ABcomp {
    tag "_${id}"
    label 'alto'
    container "dovetailg/fan-c"
    
    publishDir "${params.outDir}/AB_comp",
	mode: 'copy'
    
    input:
    tuple val(id), path(cool), val(resKB) from abcomp_mcool_ch
    	.combine(Channel.from(ABresolutions))
    
    path(genome) from abcomp_genome_ch.first()
    
    output:
    tuple val(id), path("*.bed"), path("*.ab")
    
    script:
    res = resKB.toInteger() * 1000
    """
    fanc compartments \
	-f \
	-v ${id}_eigenV_${resKB}kb.bed \
	-d ${id}_AB_${resKB}kb.bed \
	-g ${genome} \
	${cool}@${res} \
	${id}_${resKB}kb.ab
    """
    
}
