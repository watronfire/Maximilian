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
        expand( os.path.join( out_dir, "assembled/{sample}/contigs.fasta" ), out_dir=out_dir, sample=SAMPLES )

rule assemble_contigs:
    input:
        read1depleted = os.path.join( out_dir, "depleted/{sample}_R1_depleted.fastq" ),
        read2depleted = os.path.join( out_dir, "depleted/{sample}_R2_depleted.fastq" )
    output:
        os.path.join( out_dir, "assembled/{sample}/contigs.fasta" )
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
