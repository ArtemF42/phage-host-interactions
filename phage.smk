rule all:
    input:
        ...


rule seqtk:
    input:
        'reads/{sample}_{orientation}.fq.gz'
    output:
        'phage_assembly/{sample}/reads/{sample}_{orientation}.fq.gz'
    shell:
        'seqtk sample -s42 {input} 10000 | gzip > {output}'


rule spades:
    input:
        r1 = 'phage_assembly/{sample}/reads/{sample}_R1.fq.gz',
        r2 = 'phage_assembly/{sample}/reads/{sample}_R2.fq.gz'
    output:
        'phage_assembly/{sample}/assembly/contigs.fasta'
    params:
        out = 'phage_assembly/{sample}/assembly'
    threads: 256
    shell:
        'spades.py --only-assembler -t {threads} -1 {input.r1} -2 {input.r2} -o {params.out}'


rule checkv:
    input:
        'phage_assembly/{sample}/assembly/contigs.fasta'
    output:
        expand('phage_assembly/{{sample}}/checkv/{name}.tsv',
               name=['quality_summary', 'complete_genomes'])
    params:
        out = 'phage_assembly/{sample}/checkv'
    threads: 256
    conda: 'checkv'
    shell:
        'checkv end_to_end {input} {params.out} -t {threads}'


rule extraction:
    input:
        'phage_assembly/{sample}/assembly/contigs.fasta',
        expand('phage_assembly/{{sample}}/checkv/{name}.tsv',
               name=['quality_summary', 'complete_genomes'])
    output:
        'phage_assembly/{sample}/checkv/{sample}.fasta'
    run:
        import pandas as pd
        from Bio import SeqIO

        quality_summary = pd.read_csv(input[1], sep='\t', index_col='contig_id')
        quality_summary = quality_summary[quality_summary['checkv_quality'] == 'Complete']

        if len(quality_summary) == 1:
            contig_id = quality_summary.index.item()
        else:
            raise Exception('Please review the CheckV output manually')
        
        complete_genomes = pd.read_csv(input[2], sep='\t', index_col='contig_id')
        repeat_seq = complete_genomes.loc[contig_id, 'repeat_seq']

        contigs = SeqIO.index(input[0], 'fasta')

        phage_genome = contigs[contig_id]
        phage_genome.seq = phage_genome.seq.removesuffix(repeat_seq)

        with open(output[0], 'w') as file:
            SeqIO.write(phage_genome, file, 'fasta')


rule dnaapler:
    input:
        'phage_assembly/{sample}/checkv/{sample}.fasta'
    output:
        'phage_assembly/{sample}/dnaapler/{sample}_reoriented.fasta'
    params:
        out = 'phage_assembly/{sample}/dnaapler'
    threads: 256
    conda: 'pharokka'
    shell:
        'dnaapler phage -i {input} -o {params.out} -t {threads} -p {wildcards.sample} -f'


rule bowtie2:
    input:
        r1 = 'reads/{sample}_R1.fq.gz',
        r2 = 'reads/{sample}_R2.fq.gz',
        ref = 'phage_assembly/{sample}/dnaapler/{sample}_reoriented.fasta'
    output:
        'phage_assembly/{sample}/bowtie2/mapped.sam'
    params:
        idx = 'phage_assembly/{sample}/bowtie2/index/{sample}'
    threads: 256
    shell:
        'mkdir phage_assembly/{wildcards.sample}/bowtie2/index && '
        'bowtie2-build {input.ref} {params.idx} --threads {threads} && '
        'bowtie2 -x {params.idx} -1 {input.r1} -2 {input.r2} -S {output} --no-unal -p {threads}'


rule samtools:
    input:
        'phage_assembly/{sample}/bowtie2/mapped.sam'
    output:
        bam = 'phage_assembly/{sample}/bowtie2/mapped.bam',
        bai = 'phage_assembly/{sample}/bowtie2/mapped.bam.bai'
    threads: 256
    shell:
        'samtools sort -@ {threads} -o {output.bam} {input} && '
        'samtools index {output.bam}'


rule pilon:
    input:
        ref = 'phage_assembly/{sample}/dnaapler/{sample}_reoriented.fasta',
        bam = 'phage_assembly/{sample}/bowtie2/mapped.bam'
    output:
        'phage_assembly/{sample}/pilon/{sample}.fasta'
    params:
        out = 'phage_assembly/{sample}/pilon'
    conda: 'pilon'
    shell:
        'pilon --genome {input.ref} --frags {input.bam} --output {wildcards.sample} '
        '--outdir {params.out} --changes -Xmx200G'


rule pharokka:
    input:
        'phage_assembly/{sample}/pilon/{sample}.fasta'
    output:
        expand('phage_assembly/{{sample}}/pharokka/{{sample}}.{ext}',
               ext=['gbk', 'tbl'])
    params:
        db = '/home/artem42/mega/artem42/data/pharokka',
        out = 'phage_assembly/{sample}/pharokka'
    threads: 256
    conda:
        'pharokka'
    shell:
        'pharokka.py -i {input} -o {params.out} -d {params.db} -t {threads} '
        '-g prodigal -p {wildcards.sample}'
