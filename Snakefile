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
print( SAMPLES )

rule all:
    input:
        expand( os.path.join( out_dir, "trimmed/{sample}_{read}_{type}.fastq" ), out_dir=out_dir, sample=SAMPLES, read=["R1", "R2"], type=["trimmed", "unpaired"] )

rule deplete_host_reads:
    input:
        read1trimmed = os.path.join( out_dir, "trimmed/{sample}_R1_trimmed.fastq" ),
        read1unpaired = os.path.join( out_dir, "trimmed/{sample}_R1_unpaired.fastq" ),
        read2trimmed = os.path.join( out_dir, "trimmed/{sample}_R2_trimmed.fastq" ),
        read2unpaired = os.path.join( out_dir, "trimmed/{sample}_R2_unpaired.fastq" ),


rule generate_host_ref:
    input:
        host_ref = config["generate_host_ref"]["host_ref_location"]
    params:
        destination = os.path.join( config["generate_host_ref"]["host_ref_destination"], config["generate_host_ref"]["host_ref_name"] )
    output:
        host_ref_index = params.destination + ".1.bt2"
    shell:
        "bowtie2-build {input.host_ref} {output.host_ref_index}"


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
        Path( output.read1trimmed ).touch()
        Path( output.read1unpaired ).touch()
        Path( output.read2trimmed ).touch()
        Path( output.read2unpaired ).touch()
        command = "java -Xmx2g -classpath {} org.usadellab.trimmomatic.TrimmomaticPE -threads 16 {read1} {read2} {read1trimmed} {read1unpaired} {read2trimmed} {read2unpaired} {options}".format( config["adaptor_trimming"]["trimmomatic_location"],
                                                                                                                                                                                        read1=input.read1,
                                                                                                                                                                                        read2=input.read2,
                                                                                                                                                                                        read1trimmed=output.read1trimmed,
                                                                                                                                                                                        read1unpaired=output.read1unpaired,
                                                                                                                                                                                        read2trimmed=output.read2trimmed,
                                                                                                                                                                                        read2unpaired=output.read1unpaired,
                                                                                                                                                                                        options=config["adaptor_trimming"]["trimmomatic_options"])
        print( command )
        subprocess.call( command, shell=True )
