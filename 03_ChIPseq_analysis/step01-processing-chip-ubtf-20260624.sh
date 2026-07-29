
### Created on: 2026-07-16 15:58
### Created by: k2367543

cd /scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf
export out fqs fq1 fq2 fq bams bam idx

# qc

in=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/data
out=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/qc && mkdir -p ${out}
fqs=$(printf '%s\n' ${in}/* | grep -E '\.(fastq|fq)(\.gz)?$')
samples=$(basename -a ${fqs} | sed -E "s/\.(fastq|fq)(\.gz)?.*$//" | uniq)
parallel -j 3 '
    sample={}
    printf "%s\n" ${out}/* | grep -q ${sample} && exit 0
    fq=$(printf "%s\n" ${fqs} | grep ${sample})
    fastqc \
--quiet \
--outdir ${out} \
--threads 4 \
${fqs}
' ::: ${samples}

# trimming

in=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/data
out=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/trim && mkdir -p ${out}
fqs=$(printf "%s\n" ${in}/* | grep -E "\.(fastq|fq)(\.gz)?$")
samples=$(basename -a ${fqs} | sed -E "s/_R?[12][._].*$//" | uniq)
parallel -j 3 '
    sample={}
    printf "%s\n" ${out}/* | grep -q ${sample} && exit 0
    fq1=$(printf "%s\n" ${fqs} | grep ${sample} | grep -E "_R?1[._]")
    fq2=$(printf "%s\n" ${fqs} | grep ${sample} | grep -E "_R?2[._]")
    [[ -z "${fq1}" || -z "${fq2}" ]] && exit 1
    trim_galore \
--paired \
--gzip \
--cores 4 \
--output_dir ${out} \
${fq1} ${fq2}
' ::: ${samples}

# build alignment index

fasta=/scratch/prj/mmg_holland_sandbox/_genomes/mouse/GRCm38.p0/rdna/GRCm38.masked.rDNA.looped.fasta
idx=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/bowtie2-idx && mkdir -p ${idx}
bowtie2-build \
--threads 4 \
${fasta} \
${idx}/bowtie2-idx

# alignment to genome

in=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/trim
out=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/map && mkdir -p ${out}
fqs=$(printf '%s\n' ${in}/* | grep -E '\.(fastq|fq)(\.gz)?$')
samples=$(basename -a ${fqs} | sed -E "s/_R?[12][._].*$//" | uniq)
parallel -j 3 '
    sample={}
    printf "%s\n" ${out}/* | grep -q ${sample} && exit 0
    fq1=$(printf "%s\n" ${fqs} | grep ${sample} | grep -E "_R?1[._]")
    fq2=$(printf "%s\n" ${fqs} | grep ${sample} | grep -E "_R?2[._]")
    [[ -z "${fq1}" || -z "${fq2}" ]] && exit 1
    bowtie2 \
-x /scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/bowtie2-idx/bowtie2-idx \
--threads 4 \
-1 ${fq1} -2 ${fq2} \
    | samtools view \
        -b \
        --threads 4 \
    | samtools sort \
        --threads 4 \
        -o ${out}/${sample}.bam 
    samtools index ${out}/${sample}.bam
' ::: ${samples}

# deduplication

# WARNING
# If the @RG tag was not generated, the below script would generate it automatically.
# Please, modify this script manually if you require any specific @RG tag for downstream analysis.

in=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/map
out=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/dedup && mkdir -p ${out}
bams=$(printf "%s\n" ${in}/* | grep -E "\.(bam|sam)$")
samples=$(basename -a ${bams} | sed -E "s/\.(bam|sam)$//" | uniq)
parallel -j 3 '
    sample={}
    printf "%s\n" /* | grep -q  && exit 0
    bam=$(printf "%s\n" ${bams} | grep ${sample} | grep -E "\.(bam|sam)$")
    rg_tag=$(samtools view -H ${bam} | grep "^@RG" || true)
    if [[ -z "${rg_tag}" ]]; then
        picard AddOrReplaceReadGroups \
            --INPUT ${bam} \
            --OUTPUT ${out}/${sample}.rg.bam \
            --RGID ${sample} \
            --RGLB unknown \
            --RGPL ILLUMINA \
            --RGPU unknown \
            --RGSM ${sample}
        bam=${out}/${sample}.rg.bam
    fi
    sort_tag=$(samtools view -H ${bam} | grep -o "SO:[^ ]*" || true)
    if [[ "${sort_tag}" != "SO:coordinate" ]]; then
        samtools sort \
            --threads 4 \
            -o ${out}/${sample}.sort.bam \
            ${bam}
        bam=${out}/${sample}.sort.bam
    fi
    picard MarkDuplicates \
--REMOVE_DUPLICATES true \
--INPUT ${bam} \
--OUTPUT ${out}/${sample}.dedup.bam \
--METRICS_FILE ${out}/${sample}.metrics.txt 
    samtools index ${out}/${sample}.dedup.bam
    rm -f ${out}/${sample}.rg.bam
    rm -f ${out}/${sample}.sort.bam
' ::: ${samples}

# reads per region

in=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/dedup
out=/scratch/prj/mmg_holland_sandbox/Mila/_rdna_chip_ubtf/bam_BK000964_3_looped_3008 && mkdir -p ${out}
bams=$(printf '%s\n' ${in}/* | grep -E '\.(bam|sam)$')
samples=$(basename -a ${bams} | sed -E "s/\.(bam|sam).*$//" | uniq)
parallel -j 3 '
    sample={}
    printf "%s\n" ${out}/* | grep -q ${sample} && exit 0
    bam=$(printf "%s\n" ${bams} | grep ${sample})
    samtools view \
-b \
-o ${out}/${sample}.BK000964_3_looped_3008.bam \
--threads 4 \
${bam} \
BK000964.3.looped_3008 \
samtools index ${out}/${sample}.BK000964_3_looped_3008.bam
' ::: ${samples}
 