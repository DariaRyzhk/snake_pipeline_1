# rule all: основной выходной файл для завершения всех операций, указывает на финальный результат
rule all:
    input:
        "results/genes_uniq_list.txt"

rule download_gtf:
    output:
        "results/GCF_000001405.13_GRCh37_genomic.gtf"
    conda:
        "envs/wget.yaml"
    shell:
        """
        wget -O {output}.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.gtf.gz
        gunzip -c {output}.gz > {output}
        rm {output}.gz
        """

rule extract_genes_exons:
    input:
        "results/GCF_000001405.13_GRCh37_genomic.gtf"
    output:
        genes="results/genes_RefSeq_awk.gtf",
        exons="results/exons_RefSeq_awk.gtf"
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} \
            /^#/ {{print $0 > "{output.genes}"; print $0 > "{output.exons}"; next}} \
            $3 == "gene" {{print $0 > "{output.genes}"}} \
            $3 == "exon" {{print $0 > "{output.exons}"}}' {input}
        """

rule intersect_genes:
    input:
        bed="data/IAD143293_241_Designed_formatRefSeq.bed",
        gtf="results/genes_RefSeq_awk.gtf"
    output:
        "results/annotated_regions_genes_RefSeq.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.gtf} -wa -wb -nonamecheck > {output}
        """

rule extract_gene_names:
    input:
        "results/annotated_regions_genes_RefSeq.txt"
    output:
        "results/genes_list_RefSeq.txt"
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        awk '{{print $1 "\\t" $4 "\\t" substr($16, 2, length($16)-3)}}' {input} > {output}
        """

rule intersect_exons:
    input:
        bed="data/IAD143293_241_Designed_formatRefSeq.bed",
        gtf="results/exons_RefSeq_awk.gtf"
    output:
        "results/annotated_regions_genes_with_exons_RefSeq.txt"
    conda:
        "envs/bedtools.yaml"
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.gtf} -wa -wb -nonamecheck | \
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} \
            {{if ($3 == "exon") {{print $0, exon_number++}} else {{print $0, "."}}}}' > {output}
        """

rule extract_exon_info:
    input:
        "results/annotated_regions_genes_with_exons_RefSeq.txt"
    output:
        "results/genes_with_e_list_awk_RefSeq.txt"
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        awk '{{print $1 "\\t" $4 "\\t" substr($16, 2, length($16)-3) "\\t" substr($(NF-1), 2, 1)}}' {input} > {output}
        """

rule unique_gene_names:
    input:
        "results/genes_with_e_list_awk_RefSeq.txt"
    output:
        "results/genes_uniq_list.txt"
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        awk '{{print $3}}' {input} | sort -u > {output}
        """
