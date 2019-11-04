import os
from fnmatch import fnmatch
import subprocess
from pathlib import Path

in_dir = config["all"]["in_dir"]
in_ext = config["all"]["in_ext"]
out_dir = config["all"]["out_dir"]

SAMPLES = dict()
for f in os.listdir( in_dir ):
    if fnmatch( f, "*.{}".format( in_ext ) ):
        fsplit = f.split( "_" )
        if fsplit[0] not in SAMPLES:
            SAMPLES[fsplit[0]] = [f]
        else:
            SAMPLES[fsplit[0]].append( f )

rule all:
    input:
        expand( os.path.join( out_dir, "classified/{sample}.m8" ), out_dir=out_dir, sample=SAMPLES ),
        os.path.join( out_dir, "output.krona.html" ),
        [os.path.join( out_dir, "reports/{}_fastqc.html".format( os.path.basename( os.path.splitext( os.path.splitext( item )[0] )[0]     ) ) ) for sublist in SAMPLES.values() for item in sublist]

rule convert_output:
    input:
        matches = expand( os.path.join( out_dir, "classified/{sample}.m8" ), sample=SAMPLES )
    output: 
        krona_output = os.path.join( out_dir, "output.krona.html" )
    shell:
        "ktImportBLAST {config[convert_output][krona_options]}  {input.matches} -o {output.krona_output} -tax {config[convert_output][tax_location]}"

rule classify_reads:
    input:
        scaffolds_singlets = os.path.join( out_dir, "assembled/{sample}/scaffolds_singlets.fasta" )
    output:
        matches = os.path.join( out_dir, "classified/{sample}.m8" )
    shell:
        "diamond blastx -d {config[classify_reads][db_loc]} -q {input.scaffolds_singlets} -o {output.matches}" 

rule generate_singletons:
    input:
        scaffolds = os.path.join( out_dir, "assembled/{sample}/scaffolds.fasta" ),
        corrected_r1 = os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R1_depleted.00.0_0.cor.fastq.gz" ),
        corrected_r2 = os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R2_depleted.00.0_0.cor.fastq.gz" ),
        corrected_unpaired = os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R_unpaired.00.0_0.cor.fastq.gz" )
    output:
        scaffolds_singlets = os.path.join( out_dir, "assembled/{sample}/scaffolds_singlets.fasta" ),
        scaffold_aligned = temp( os.path.join( out_dir, "assembled/{sample}/scaffold_aligned.bam" ) ),
        scaffold_unmapped = temp( os.path.join( out_dir, "assembled/{sample}/scaffold_unmapped.bam" ) ),
        singlets_fasta = temp( os.path.join( out_dir, "assembled/{sample}/singlets.fasta" ) )
    run:
        depleted_index = os.path.join( os.path.dirname( input.scaffolds ), "scaffolds" )
        build_command = "bowtie2-build {} {}".format( input.scaffolds, depleted_index )        
        
        map_command = "bowtie2 -x {} -1 {} -2 {} -U {} | samtools view -bS - > {}".format( depleted_index, input.corrected_r1, input.corrected_r2, input.corrected_unpaired, output.scaffold_aligned )

        filter_command = "samtools view -b -f 12 -F 256 {} > {}".format( output.scaffold_aligned, output.scaffold_unmapped )
        
        fasta_command = "samtools fasta {} > {}".format( output.scaffold_unmapped, output.singlets_fasta )
        
        merge_command = "cat {} {} > {}".format( input.scaffolds, output.singlets_fasta, output.scaffolds_singlets )       

        linked_command = " && ".join( ["module load samtools",
                                       "module load bowtie2",
                                       build_command,
                                       map_command,
                                       filter_command,
                                       fasta_command,
                                       merge_command] )

        subprocess.call( linked_command, shell=True )
rule assemble_contigs:
    input:
        read1depleted = os.path.join( out_dir, "depleted/{sample}_R1_depleted.fastq" ),
        read2depleted = os.path.join( out_dir, "depleted/{sample}_R2_depleted.fastq" )
    output:
        os.path.join( out_dir, "assembled/{sample}/scaffolds.fasta" ),
        os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R1_depleted.00.0_0.cor.fastq.gz" ),
        os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R2_depleted.00.0_0.cor.fastq.gz" ),
        os.path.join( out_dir, "assembled/{sample}/corrected/{sample}_R_unpaired.00.0_0.cor.fastq.gz" )
    run:
        command = "spades.py --meta -k 21,33,55,77 -1 {} -2 {} -o {}".format( input.read1depleted, input.read2depleted, os.path.join( out_dir, "assembled", wildcards.sample ) )
        subprocess.call( command, shell=True ) 

rule deplete_host_reads:
    input:
        read1trimmed = os.path.join( out_dir, "trimmed/{sample}_R1_trimmed.fastq" ),
        read1unpaired = os.path.join( out_dir, "trimmed/{sample}_R1_unpaired.fastq" ),
        read2trimmed = os.path.join( out_dir, "trimmed/{sample}_R2_trimmed.fastq" ),
        read2unpaired = os.path.join( out_dir, "trimmed/{sample}_R2_unpaired.fastq" ),
        host_ref = os.path.join( config["generate_host_ref"]["host_ref_destination"], "{}.1.bt2".format( config["generate_host_ref"]["host_ref_name"] ) )
    params:
        host_ref = os.path.join( config["generate_host_ref"]["host_ref_destination"], config["generate_host_ref"]["host_ref_name"] )
    output:
        host_aligned = temp( os.path.join( out_dir, "depleted/{sample}_aligned.bam" ) ),
        unmapped = temp( os.path.join( out_dir, "depleted/{sample}_unmapped.bam" ) ),
        sorteed = temp( os.path.join( out_dir, "depleted/{sample}_sorted.bam" ) ),
        read1depleted = os.path.join( out_dir, "depleted/{sample}_R1_depleted.fastq" ),
        read2depleted = os.path.join( out_dir, "depleted/{sample}_R2_depleted.fastq" )
    shell:
        "module load bowtie2 &&"
        "module load bedtools &&"
        "bowtie2 -x {params.host_ref} {config[deplete_host_reads][bowtie_options]} -1 {input.read1trimmed} -2 {input.read2trimmed} -U {input.read1unpaired},{input.read2unpaired} | samtools view -bS - > {output.host_aligned} &&"
        "samtools view -b -f 12 -F 256 {output.host_aligned} > {output.unmapped} &&"
        "samtools sort -n {output.unmapped} -o {output.sorteed} &&"
        "bedtools bamtofastq -i {output.sorteed} -fq {output.read2depleted} -fq2 {output.read1depleted}"

rule generate_host_ref:
    input:
        host_ref = config["generate_host_ref"]["host_ref_location"]
    params:
        destination = os.path.join( config["generate_host_ref"]["host_ref_destination"], config["generate_host_ref"]["host_ref_name"] )
    output:
        host_ref_index = os.path.join( config["generate_host_ref"]["host_ref_destination"], "{}.1.bt2".format( config["generate_host_ref"]["host_ref_name"] ) )
    shell:
        "bowtie2-build {input.host_ref} {params.destination}"


rule trim_adaptors:
    input:
        read1 = lambda wildcards: os.path.join( "{in_dir}".format( in_dir = in_dir ), SAMPLES[wildcards.sample][0] ),
        read2 = lambda wildcards: os.path.join( "{in_dir}".format( in_dir = in_dir ), SAMPLES[wildcards.sample][1] )
    output:
        read1trimmed = os.path.join( out_dir, "trimmed/{sample}_R1_trimmed.fastq" ),
        read1unpaired = os.path.join( out_dir, "trimmed/{sample}_R1_unpaired.fastq" ),
        read2trimmed = os.path.join( out_dir, "trimmed/{sample}_R2_trimmed.fastq" ),
        read2unpaired = os.path.join( out_dir, "trimmed/{sample}_R2_unpaired.fastq" )
    run:
        command = "java -Xmx2g -classpath {} org.usadellab.trimmomatic.TrimmomaticPE -threads 16 {read1} {read2} {read1trimmed} {read1unpaired} {read2trimmed} {read2unpaired} {options}".format( config["adaptor_trimming"]["trimmomatic_location"],
                                                                                                                                                                                        read1=input.read1,
                                                                                                                                                                                        read2=input.read2,
                                                                                                                                                                                        read1trimmed=output.read1trimmed,
                                                                                                                                                                                        read1unpaired=output.read1unpaired,
                                                                                                                                                                                        read2trimmed=output.read2trimmed,
                                                                                                                                                                                        read2unpaired=output.read2unpaired,
                                                                                                                                                                                        options=config["adaptor_trimming"]["trimmomatic_options"] )
        subprocess.call( command, shell=True )
        Path( output.read1trimmed ).touch()
        Path( output.read1unpaired ).touch()
        Path( output.read2trimmed ).touch()
        Path( output.read2unpaired ).touch()

rule quality_control:
    input:
        reads = [os.path.join( "{in_dir}".format( in_dir = in_dir ), item ) for sublist in SAMPLES.values() for item in sublist]
    params:
        reports_location = os.path.join( out_dir, "reports" )
    output:
        reports = [os.path.join( out_dir, "reports/{}_fastqc.html".format( os.path.basename( os.path.splitext( os.path.splitext( item )[0] )[0] ) ) ) for sublist in SAMPLES.values() for item in sublist]
    shell:
        "module load fastqc &&"
        "fastqc -o {params.reports_location} --nogroup {input.reads}"
