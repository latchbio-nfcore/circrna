from typing import List, Optional

from latch.resources.launch_plan import LaunchPlan
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile

from wf.entrypoint import Reference, SampleSheet, initialize, nextflow_runtime


@workflow(metadata._nextflow_metadata)
def nf_nf_core_circrna(
    run_name: str,
    input: List[SampleSheet],
    phenotype: Optional[LatchFile],
    annotation: Optional[LatchFile],
    email: Optional[str],
    multiqc_title: Optional[str],
    mirna_expression: Optional[LatchFile],
    seq_center: Optional[str],
    genome_source: str,
    fasta: Optional[LatchFile],
    gtf: Optional[LatchFile],
    mature: Optional[LatchFile],
    bowtie: Optional[LatchDir],
    bowtie2: Optional[LatchDir],
    bwa: Optional[LatchDir],
    hisat2: Optional[LatchDir],
    segemehl: Optional[LatchFile],
    star: Optional[LatchDir],
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    multiqc_methods_description: Optional[str],
    tools: str = "circexplorer2",
    bsj_reads: Optional[int] = 1,
    max_shift: Optional[int] = 1,
    min_tools: Optional[int] = 1,
    min_samples: Optional[int] = 1,
    exon_boundary: Optional[int] = 0,
    quantification_tools: Optional[str] = "ciriquant,psirc",
    bootstrap_samples: Optional[int] = 30,
    mirna_min_sample_percentage: Optional[float] = 0.2,
    mirna_min_reads: Optional[int] = 5,
    mirna_correlation: Optional[str] = "pearson",
    mirna_tools: Optional[str] = "miranda,targetscan",
    mirna_min_tools: Optional[int] = 1,
    sjdboverhang: Optional[int] = 100,
    chimJunctionOverhangMin: Optional[int] = 10,
    alignSJDBoverhangMin: Optional[int] = 10,
    limitSjdbInsertNsj: Optional[int] = 1000000,
    chimSegmentMin: Optional[int] = 10,
    seglen: Optional[int] = 25,
    min_intron: Optional[int] = 20,
    max_intron: Optional[int] = 1000000,
    min_map_len: Optional[int] = 40,
    min_fusion_distance: Optional[int] = 200,
    save_unaligned: bool = False,
    save_reference: bool = True,
    hisat2_build_memory: Optional[str] = "200.GB",
    skip_trimming: bool = False,
    save_trimmed: bool = False,
    skip_fastqc: bool = False,
    min_trimmed_reads: Optional[int] = 10000,
    save_intermediates: bool = False,
    genome: Reference = Reference.GRCh38,
    outdir: LatchOutputDir = LatchOutputDir("latch:///CircRNA"),
) -> None:
    """
    nf-core/circrna is a bioinformatics pipeline to analyse total RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming, back-splice junction (BSJ) detection, annotation, quantification and miRNA target prediction of circular RNAs.

    <html>
    <p align="center">
    <img src="https://user-images.githubusercontent.com/31255434/182289305-4cc620e3-86ae-480f-9b61-6ca83283caa5.jpg" alt="Latch Verified" width="100">
    </p>

    <p align="center">
    <strong>
    Latch Verified
    </strong>
    </p>

    <p align="center">

    ## Introduction

    [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

    **nf-core/circrna** is a bioinformatics pipeline to analyse total RNA sequencing data obtained from organisms with a reference genome and annotation. It takes a samplesheet and FASTQ files as input, performs quality control (QC), trimming, back-splice junction (BSJ) detection, annotation, quantification and miRNA target prediction of circular RNAs.

    The pipeline is still under development, but the BSJ detection and quantification steps are already implemented and functional. The following features are planned to be implemented soon:

    - Isoform-level circRNA detection and quantification
    - circRNA-miRNA interaction analysis using [SPONGE](https://doi.org/10.1093/bioinformatics/btz314) and [spongEffects](https://doi.org/10.1093/bioinformatics/btad276)
    - Improved downstream analyses

    If you want to contribute, feel free to create an issue or pull request on the [GitHub repository](https://github.com/nf-core/circrna) or join the [Slack channel](https://nf-co.re/join/slack).

    ## Pipeline summary

    - Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
    - Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
    - BSJ detection
    - [`CIRIquant`](https://github.com/Kevinzjy/CIRIquant)
    - [`STAR 2-Pass mode`](https://github.com/alexdobin/STAR)
        - [`CIRCexplorer2`](https://circexplorer2.readthedocs.io/en/latest/)
        - [`circRNA finder`](https://github.com/orzechoj/circRNA_finder)
        - [`DCC`](https://github.com/dieterich-lab/DCC)
    - [`find circ`](https://github.com/marvin-jens/find_circ)
    - [`MapSplice`](http://www.netlab.uky.edu/p/bioinfo/MapSplice2)
    - [`Segemehl`](https://www.bioinf.uni-leipzig.de/Software/segemehl/)
    - circRNA annotation
    - Based on a GTF file
    - Based on database files (if provided)
    - Extract circRNA sequences and build circular transcriptome
    - Merge circular transcriptome with linear transcriptome derived from provided GTF
    - Quantification of combined circular and linear transcriptome
    - [`psirc-quant`](https://github.com/Christina-hshi/psirc)
    - miRNA binding affinity analysis (only if the `mature` parameter is provided)
    - Normalizes miRNA expression (only if the `mirna_expression` parameter is provided)
    - Binding site prediction
        - [`miRanda`](http://cbio.mskcc.org/miRNA2003/miranda.html)
        - [`TargetScan`](http://www.targetscan.org/cgi-bin/targetscan/data_download.vert72.cgi)
    - Perform majority vote on binding sites
    - Compute correlations between miRNA and transcript expression levels (only if the `mirna_expression` parameter is provided)
    - Statistical tests (only if the `phenotype` parameter is provided)
    - [`CircTest`](https://github.com/dieterich-lab/CircTest)
    - MultiQC report [`MultiQC`](http://multiqc.info/)

    ## Usage

    First, prepare a samplesheet with your input data that looks as follows:

    ```
    sample,fastq_1,fastq_2
    CONTROL,CONTROL_R1.fastq.gz,CONTROL_R2.fastq.gz
    TREATMENT,TREATMENT_R1.fastq.gz,TREATMENT_R2.fastq.gz
    ```

    Each row represents a fastq file (single-end) or a pair of fastq files (paired end).
    For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/circrna/usage) and the [parameter documentation](https://nf-co.re/circrna/parameters).

    ## Pipeline output

    To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/circrna/results) tab on the nf-core website pipeline page.
    For more details about the output files and reports, please refer to the
    [output documentation](https://nf-co.re/circrna/output).

    ## Credits

    nf-core/circrna was originally written by [Barry Digby](https://github.com/BarryDigby).
    It was later refactored, extended and improved by [Nico Trummer](https://github.com/nictru).

    We thank the following people for their extensive assistance in the development of this pipeline (in alphabetical order):

    - [Alexander Peltzer](https://github.com/apeltzer)
    - [Ben Whittle](https://github.com/bj-w)
    - [Kevin Menden](https://github.com/KevinMenden)
    - [Malte Weyrich](https://github.com/mweyrich28)
    - [Marieke Vromman](https://github.com/MariekeVromman)
    - [Maxime Garcia](https://github.com/maxulysse)
    - [Phil Ewels](https://github.com/ewels)

    ## Acknowledgements

    ![SFI](./docs/images/Genomics-Data-Science-original.png)

    ## Contributions and Support

    If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

    For further information or help, don't hesitate to get in touch on the [Slack `#circrna` channel](https://nfcore.slack.com/channels/circrna) (you can join with [this invite](https://nf-co.re/join/slack)).

    ## Citations

    <!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
    <!-- If you use nf-core/circrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

    > **nf-core/circrna: a portable workflow for the quantification, miRNA target prediction and differential expression analysis of circular RNAs.**
    >
    > Barry Digby, Stephen P. Finn, & Pilib Ó Broin
    >
    > [BMC Bioinformatics 24, 27 (2023)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05125-8)
    > doi: [10.1186/s12859-022-05125-8](https://doi.org/10.1186/s12859-022-05125-8)

    An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file on the Github repository.

    You can cite the `nf-core` publication as follows:

    > **The nf-core framework for community-curated bioinformatics pipelines.**
    >
    > Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
    >
    > _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
        run_name=run_name,
        pvc_name=pvc_name,
        input=input,
        outdir=outdir,
        phenotype=phenotype,
        annotation=annotation,
        email=email,
        multiqc_title=multiqc_title,
        tools=tools,
        bsj_reads=bsj_reads,
        max_shift=max_shift,
        min_tools=min_tools,
        min_samples=min_samples,
        exon_boundary=exon_boundary,
        quantification_tools=quantification_tools,
        bootstrap_samples=bootstrap_samples,
        mirna_expression=mirna_expression,
        mirna_min_sample_percentage=mirna_min_sample_percentage,
        mirna_min_reads=mirna_min_reads,
        mirna_correlation=mirna_correlation,
        mirna_tools=mirna_tools,
        mirna_min_tools=mirna_min_tools,
        sjdboverhang=sjdboverhang,
        chimJunctionOverhangMin=chimJunctionOverhangMin,
        alignSJDBoverhangMin=alignSJDBoverhangMin,
        limitSjdbInsertNsj=limitSjdbInsertNsj,
        chimSegmentMin=chimSegmentMin,
        seglen=seglen,
        min_intron=min_intron,
        max_intron=max_intron,
        min_map_len=min_map_len,
        min_fusion_distance=min_fusion_distance,
        seq_center=seq_center,
        save_unaligned=save_unaligned,
        save_reference=save_reference,
        genome_source=genome_source,
        genome=genome,
        fasta=fasta,
        gtf=gtf,
        mature=mature,
        bowtie=bowtie,
        bowtie2=bowtie2,
        bwa=bwa,
        hisat2=hisat2,
        hisat2_build_memory=hisat2_build_memory,
        segemehl=segemehl,
        star=star,
        skip_trimming=skip_trimming,
        save_trimmed=save_trimmed,
        skip_fastqc=skip_fastqc,
        clip_r1=clip_r1,
        clip_r2=clip_r2,
        three_prime_clip_r1=three_prime_clip_r1,
        three_prime_clip_r2=three_prime_clip_r2,
        trim_nextseq=trim_nextseq,
        min_trimmed_reads=min_trimmed_reads,
        save_intermediates=save_intermediates,
        multiqc_methods_description=multiqc_methods_description,
    )


LaunchPlan(
    nf_nf_core_circrna,
    "Test Data",
    {
        "input": [
            SampleSheet(
                sample="cel_N2_1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep1_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep1_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="cel_N2_2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep2_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep2_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="cel_N2_3",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep3_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/N2_rep3_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="fust1_1",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep1_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep1_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="fust1_2",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep2_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep2_2.fastq.gz"
                ),
                strandedness="auto",
            ),
            SampleSheet(
                sample="fust1_3",
                fastq_1=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep3_1.fastq.gz"
                ),
                fastq_2=LatchFile(
                    "s3://latch-public/nf-core/circrna/test_data/fust1_rep3_2.fastq.gz"
                ),
                strandedness="auto",
            ),
        ],
        "genome_source": "custom",
        "fasta": LatchFile("s3://latch-public/nf-core/circrna/test_data/chrI.fa"),
        "gtf": LatchFile("s3://latch-public/nf-core/circrna/test_data/chrI.gtf"),
        "run_name": "Test_Run",
        "tools": "circexplorer2",
        "skip_trimming": False,
        "bsj_reads": 2,
    },
)
