rule all:
    input:
        "results/genes_uniq_list.txt"  # Это конечный файл, который является результатом выполнения пайплайна

rule download_gtf:
    output:
        "results/GCF_000001405.13_GRCh37_genomic.gtf"  # Скачиваем GTF-файл с аннотацией генома
    conda:
        "envs/wget.yaml"  # Используем conda окружение для скачивания
    shell:
        """
        wget -O {output}.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/GCF_000001405.13_GRCh37_genomic.gtf.gz
        gunzip -c {output}.gz > {output}
        rm {output}.gz
        """  # Скачивание GTF-файла, разархивирование и удаление архива

rule extract_genes_exons:
    input:
        "results/GCF_000001405.13_GRCh37_genomic.gtf"  # Исходный GTF-файл с аннотацией
    output:
        genes="results/genes_RefSeq_awk.gtf",  # Файл с аннотацией генов
        exons="results/exons_RefSeq_awk.gtf"  # Файл с аннотацией экзонов
    conda:
        "envs/coreutils.yaml"  # Используется conda окружение для работы с утилитами
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} \
            /^#/ {{print $0 > "{output.genes}"; print $0 > "{output.exons}"; next}} \
            $3 == "gene" {{print $0 > "{output.genes}"}} \
            $3 == "exon" {{print $0 > "{output.exons}"}}' {input}
        """  # Фильтрация GTF-файла для разделения генов и экзонов в отдельные файлы

rule intersect_genes:
    input:
        bed="data/IAD143293_241_Designed_formatRefSeq.bed",  # Файл с диапазонами для аннотации генов
        gtf="results/genes_RefSeq_awk.gtf"  # Файл с аннотацией генов
    output:
        "results/annotated_regions_genes_RefSeq.txt"  # Аннотированные регионы с генами
    conda:
        "envs/bedtools.yaml"  # Используется conda окружение для работы с bedtools
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.gtf} -wa -wb -nonamecheck > {output}
        """  # Пересечение диапазонов из .bed с аннотацией генов для получения аннотированных регионов

rule extract_gene_names:
    input:
        "results/annotated_regions_genes_RefSeq.txt"  # Аннотированные регионы с генами
    output:
        "results/genes_list_RefSeq.txt"  # Список генов из аннотированных данных
    conda:
        "envs/coreutils.yaml"  # Используется conda окружение для работы с утилитами
    shell:
        """
        awk '{{print $1 "\\t" $4 "\\t" substr($16, 2, length($16)-3)}}' {input} > {output}
        """  # Извлечение информации о генах из аннотированных данных: хромосома, имя гена и идентификатор

rule intersect_exons:
    input:
        bed="data/IAD143293_241_Designed_formatRefSeq.bed",  # Файл с диапазонами для аннотации экзонов
        gtf="results/exons_RefSeq_awk.gtf"  # Файл с аннотацией экзонов
    output:
        "results/annotated_regions_genes_with_exons_RefSeq.txt"  # Аннотированные регионы с экзонами
    conda:
        "envs/bedtools.yaml"  # Используется conda окружение для работы с bedtools
    shell:
        """
        bedtools intersect -a {input.bed} -b {input.gtf} -wa -wb -nonamecheck | \
        awk 'BEGIN {{FS="\\t"; OFS="\\t"}} \
            {{if ($3 == "exon") {{print $0, exon_number++}} else {{print $0, "."}}}}' > {output}
        """  # Пересечение диапазонов с экзонами и добавление номера экзона в результат

rule extract_exon_info:
    input:
        "results/annotated_regions_genes_with_exons_RefSeq.txt"  # Аннотированные регионы с экзонами
    output:
        "results/genes_with_e_list_awk_RefSeq.txt"  # Информация о генах с экзонами
    conda:
        "envs/coreutils.yaml"  # Используется conda окружение для работы с утилитами
    shell:
        """
        awk '{{print $1 "\\t" $4 "\\t" substr($16, 2, length($16)-3) "\\t" substr($(NF-1), 2, 1)}}' {input} > {output}
        """  # Извлечение информации о генах и экзонах: хромосома, имя гена, экзон, номер экзона

rule unique_gene_names:
    input:
        "results/genes_with_e_list_awk_RefSeq.txt"  # Информация о генах с экзонами
    output:
        "results/genes_uniq_list.txt"  # Список уникальных генов
    conda:
        "envs/coreutils.yaml"  # Используется conda окружение для работы с утилитами
    shell:
        """
        awk '{{print $3}}' {input} | sort -u > {output}
        """  # Извлечение уникальных имен генов из списка и сортировка
