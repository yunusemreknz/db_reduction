SAMPLES = [""]  # Base names of FASTA files WITHOUT ".fasta" extension

rule all:
    input:
        expand("reducedfasta/{sample}_filtered.fasta", sample=SAMPLES)

rule insilicodigestion:
    input:
        fasta="data/{sample}.fasta"
    output:
        digested="data/{sample}_trypsin_digestion.txt"
    conda:
        "deepmspeptide2"
    shell:
        "python scripts/insilicodigestion.py {input.fasta} {output.digested}"

rule dmsp:
    input:
        infile="data/{sample}_trypsin_digestion.txt",
        model="model/model_2_1D.h5"
    output:
        predictions="data/{sample}_trypsin_digestion.Predictions.txt"
    conda:
        "deepmspeptide2"
    shell:
        "python scripts/dmsp.py {input.infile} {output.predictions}"

rule createfasta:
    input:
        predictions="data/{sample}_trypsin_digestion.Predictions.txt"
    output:
        fasta="reducedfasta/{sample}_filtered.fasta"
    conda:
        "deepmspeptide2"
    shell:
        "python scripts/createfasta.py {input.predictions} {output.fasta}"
