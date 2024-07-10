import pdb
import pandas as pd
import os
import re
import glob

with open('sampleNames.txt', 'r') as f:
    sample_names = [line.strip() for line in f]


# Default rule
rule all:
    input:
        'DNX_data_count_anno.tsv'
        #expand('/home/mwang/dnx_projects/EXP23003664_aged_HDF/data_bedGraph/{sample}.bed', sample = sample_names)
        #expand('peakCountBySample_broadPeaks_featureCounts/{sample}.counts.txt', sample = sample_names)
        #expand('/home/mwang/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks_gtf/{sample}.gtf', sample = sample_names)
        #expand('/home/mwang/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks/{sample}.plus/MACS_RNA', sample = sample_names)
        #expand('/home/mwang/kidney/multiome/unannotatedRNA/STAR_out_allPeaksMerge_100/{sample}Solo.out/Gene/filtered/matrix.mtx', sample >

# rule to process bam files into bedGRaph
rule makeBedGraph_pos:
    input: 'data/{sample}.bam'
    output: temp('data_bedGraph/{sample}.plus.bedGraph')
    shell: """ ~/tools/bedtools genomecov -ibam {input} -bg -strand + -pc > {output} """

rule makeBedGraph_neg:
    input: 'data/{sample}.bam'
    output: temp('data_bedGraph/{sample}.neg.bedGraph')
    shell: """ ~/tools/bedtools genomecov -ibam {input} -bg -strand - > {output} """

rule combineBothDir:
    input: pos = 'data_bedGraph/{sample}.plus.bedGraph',
           neg = 'data_bedGraph/{sample}.neg.bedGraph'
    output: pos = temp('data_bedGraph/{sample}.plus.bed'),
            neg = temp('data_bedGraph/{sample}.neg.bed')
    shell: """Rscript makeBedCombine.R {input.pos} {input.neg} {output.pos} {output.neg}"""

# Rule to process each sample
rule call_peaks:
    input: pos = 'data_bedGraph/{sample}.plus.bed',
           neg = 'data_bedGraph/{sample}.neg.bed'
    output:
        pos = directory('peakCountBySample_broadPeaks/{sample}.plus/MACS_RNA'),
        neg = directory('peakCountBySample_broadPeaks/{sample}.neg/MACS_RNA')
    shell: """
        macs3 callpeak --broad -t {input.pos} --qvalue 0.01 --format BED --outdir {output.pos} -g hs --nomodel;
        macs3 callpeak --broad -t {input.neg} --qvalue 0.01 --format BED --outdir {output.neg} -g hs --nomodel
	"""

rule generatePeaksGtf:
    input:
        pos = 'peakCountBySample_broadPeaks/{sample}.plus/MACS_RNA',
        neg = 'peakCountBySample_broadPeaks/{sample}.neg/MACS_RNA'
    output:
        'peakCountBySample_broadPeaks_gtf/{sample}.gtf'
    shell: """
        Rscript processDNX.R {input.neg}/NA_peaks.xls {input.pos}/NA_peaks.xls {output}
        """

rule countPeaks:
    input:
        anno = 'peakCountBySample_broadPeaks_gtf/{sample}.gtf',
        bam = 'data/{sample}.bam'
    output:
        'peakCountBySample_broadPeaks_featureCounts/{sample}.counts.txt'
    threads: 4
    shell: """
        featureCounts -M --fraction -T {threads} -p --countReadPairs -g gene_id -t sequence_feature -s 1 -a {input.anno} -o {output} {input.bam} # highest percentage s -1
        """

rule makeDNXMat:
    input:
        sampleCounts = expand('peakCountBySample_broadPeaks_featureCounts/{sample}.counts.txt',sample=sample_names)
    output: 'DNX_data_count_anno.tsv'
    params: 'peakCountBySample_broadPeaks_featureCounts'
    threads: 4
    shell: """
        Rscript analyzeDNX_pipeline.R {params} {output}
        """
