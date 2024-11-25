from dataclasses import dataclass
from enum import Enum
from typing import List, Optional

from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch.types.metadata import (
    Fork,
    ForkBranch,
    LatchRule,
    NextflowParameter,
    Params,
    Section,
    Spoiler,
    Text,
)


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    strandedness: Optional[str]


flow = [
    Section(
        "Inputs",
        Params("input"),
        Text(
            "Sample identifier and FASTQ files should not contain spaces in file names or full directory locations."
        ),
        Text(
            "Strandedness can be set to 'auto', 'reverse', 'forward'. If left untoggled, it will default to 'auto'."
        ),
        Params(
            "phenotype",
            "annotation",
        ),
    ),
    Section(
        "Reference Genome",
        Fork(
            "genome_source",
            "",
            igenome=ForkBranch(
                "iGenome",
                Params(
                    "genome",
                ),
            ),
            custom=ForkBranch(
                "Custom Reference Genome",
                Params(
                    "fasta",
                    "gtf",
                ),
                Spoiler(
                    "Additional Reference Options",
                    Text("Paths to pre-built indices for various aligners."),
                    Params(
                        "bowtie",
                        "bowtie2",
                        "bwa",
                        "hisat2",
                        "hisat2_build_memory",
                        "segemehl",
                        "star",
                    ),
                ),
            ),
        ),
    ),
    Section(
        "Output Directory",
        Params("run_name"),
        Text("Parent directory for outputs"),
        Params("outdir"),
    ),
    Spoiler(
        "CircRNA Configuration",
        Text("Parameters controlling circRNA detection and filtering criteria."),
        Params(
            "tools",
            "bsj_reads",
            "max_shift",
            "min_tools",
            "min_samples",
            "exon_boundary",
        ),
    ),
    Spoiler(
        "Quantification",
        Text("Tools and parameters for circRNA quantification."),
        Params(
            "quantification_tools",
            "bootstrap_samples",
        ),
    ),
    Spoiler(
        "miRNA Analysis",
        Text("Parameters for miRNA expression analysis and binding site prediction."),
        Params(
            "mirna_expression",
            "mirna_min_sample_percentage",
            "mirna_min_reads",
            "mirna_correlation",
            "mirna_tools",
            "mirna_min_tools",
            "mature",
        ),
    ),
    Spoiler(
        "Alignment Options",
        Text("Parameters controlling read alignment and splice junction detection."),
        Params(
            "sjdboverhang",
            "chimJunctionOverhangMin",
            "alignSJDBoverhangMin",
            "limitSjdbInsertNsj",
            "chimSegmentMin",
            "seglen",
            "min_intron",
            "max_intron",
            "min_map_len",
            "min_fusion_distance",
            "seq_center",
            "save_unaligned",
        ),
    ),
    Spoiler(
        "Read Processing",
        Text("Options for read trimming and quality control."),
        Params(
            "skip_trimming",
            "save_trimmed",
            "skip_fastqc",
            "clip_r1",
            "clip_r2",
            "three_prime_clip_r1",
            "three_prime_clip_r2",
            "trim_nextseq",
            "min_trimmed_reads",
        ),
    ),
    Spoiler(
        "Additional Options",
        Text("Optional parameters for additional analyses and output control."),
        Params(
            # "email",
            "multiqc_title",
            "save_intermediates",
            "multiqc_methods_description",
        ),
    ),
]


generated_parameters = {
    "run_name": NextflowParameter(
        display_name="Run Name",
        description="Name of run",
        batch_table_column=True,
        rules=[
            LatchRule(
                regex=r"^[a-zA-Z0-9_-]+$",
                message="Run name must contain only letters, digits, underscores, and dashes. No spaces are allowed.",
            )
        ],
    ),
    "input": NextflowParameter(
        type=List[SampleSheet],
        display_name="Input Samples",
        description="Path to comma-separated file containing information about the samples in the experiment.",
        samplesheet=True,
        samplesheet_constructor="csv",
    ),
    "outdir": NextflowParameter(
        type=LatchOutputDir,
        display_name="Output Directory",
        description="The output directory where the results will be saved.",
    ),
    "phenotype": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Phenotype Data",
        description="Phenotype CSV file specifying the experimental design. If provided, the pipeline will run CIRCTEST.",
    ),
    "annotation": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Annotation BED Files",
        description="Path to a CSV file containing BED files that should be used for annotation.",
    ),
    "email": NextflowParameter(
        type=Optional[str],
        display_name="Notification Email",
        description="Email address for completion summary.",
    ),
    "multiqc_title": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Report Title",
        description="MultiQC report title. Printed as page header, used for filename if not otherwise specified.",
    ),
    "tools": NextflowParameter(
        type=str,
        display_name="CircRNA Detection Tools",
        description="Comma separated list of circRNA quantification tools to use: ciriquant, circexplorer2, find_circ, circrna_finder, mapsplice, dcc, segemehl",
        default="circexplorer2",
    ),
    "bsj_reads": NextflowParameter(
        type=Optional[int],
        display_name="Min Back-Splice Junction Reads",
        description="Minimum number of reads spanning circRNA back-splice junction required.",
        default=1,
    ),
    "max_shift": NextflowParameter(
        type=Optional[int],
        display_name="Max BSJ Shift",
        description="If both start and end of a pair of BSJs are within max_shift bp, they are considered as the same BSJ.",
        default=1,
    ),
    "min_tools": NextflowParameter(
        type=Optional[int],
        display_name="Min Required Tools",
        description="Minimum number of tools that must detect a circRNA for it to be included in output.",
        default=1,
    ),
    "min_samples": NextflowParameter(
        type=Optional[int],
        display_name="Min Required Samples",
        description="Minimum number of samples a circRNA must be detected in.",
        default=1,
    ),
    "exon_boundary": NextflowParameter(
        type=Optional[int],
        display_name="Exon Boundary Distance",
        description="Distance threshold for determining if a candidate is a circRNA or EI-circRNA.",
    ),
    "quantification_tools": NextflowParameter(
        type=Optional[str],
        display_name="Quantification Tools",
        description="Comma separated list of circRNA quantification tools: ciriquant, psirc",
        default="ciriquant,psirc",
    ),
    "bootstrap_samples": NextflowParameter(
        type=Optional[int],
        display_name="Bootstrap Samples",
        description="Number of bootstrap samples to use during psirc quantification.",
        default=30,
    ),
    "mirna_expression": NextflowParameter(
        type=Optional[LatchFile],
        display_name="miRNA Expression Data",
        description="Tab-separated file with miRNA expression counts from smrnaseq pipeline.",
    ),
    "mirna_min_sample_percentage": NextflowParameter(
        type=Optional[float],
        display_name="Min Sample Percentage (miRNA)",
        description="Minimum percentage of samples a miRNA must be expressed in.",
        default=0.2,
    ),
    "mirna_min_reads": NextflowParameter(
        type=Optional[int],
        display_name="Min miRNA Reads",
        description="Minimum number of reads required for a miRNA.",
        default=5,
    ),
    "mirna_correlation": NextflowParameter(
        type=Optional[str],
        display_name="miRNA Correlation Method",
        description="Type of correlation for miRNA-transcript expression analysis: pearson or spearman",
        default="pearson",
    ),
    "mirna_tools": NextflowParameter(
        type=Optional[str],
        display_name="miRNA Prediction Tools",
        description="Comma separated list of miRNA binding site prediction tools: miranda, targetscan",
        default="miranda,targetscan",
    ),
    "mirna_min_tools": NextflowParameter(
        type=Optional[int],
        display_name="Min miRNA Tools",
        description="Minimum number of tools required to predict a miRNA interaction.",
        default=1,
    ),
    "sjdboverhang": NextflowParameter(
        type=Optional[int],
        display_name="STAR Junction Overhang",
        description="Number of bases to use from donor/acceptor sides of junctions for STAR indexing.",
        default=100,
    ),
    "chimJunctionOverhangMin": NextflowParameter(
        type=Optional[int],
        display_name="Min Chimeric Junction Overhang",
        description="Minimum overhang for a chimeric junction.",
        default=10,
    ),
    "alignSJDBoverhangMin": NextflowParameter(
        type=Optional[int],
        display_name="Min Junction Database Overhang",
        description="Minimum overhang for annotated junctions.",
        default=10,
    ),
    "limitSjdbInsertNsj": NextflowParameter(
        type=Optional[int],
        display_name="Max Junction Insertions",
        description="Maximum number of junctions to insert into genome during mapping.",
        default=1000000,
    ),
    "chimSegmentMin": NextflowParameter(
        type=Optional[int],
        display_name="Min Chimeric Segment Length",
        description="Minimum length of chimeric segment. Must be positive for circular detection.",
        default=10,
    ),
    "seglen": NextflowParameter(
        type=Optional[int],
        display_name="Segment Length",
        description="Length of segments for alignment.",
        default=25,
    ),
    "min_intron": NextflowParameter(
        type=Optional[int],
        display_name="Min Intron Length",
        description="Minimum intron length for splicing.",
        default=20,
    ),
    "max_intron": NextflowParameter(
        type=Optional[int],
        display_name="Max Intron Length",
        description="Maximum intron length for splicing.",
        default=1000000,
    ),
    "min_map_len": NextflowParameter(
        type=Optional[int],
        display_name="Min Mapping Length",
        description="Minimum alignment length required.",
        default=40,
    ),
    "min_fusion_distance": NextflowParameter(
        type=Optional[int],
        display_name="Min Fusion Distance",
        description="Minimum distance between gapped segments for fusion candidates.",
        default=200,
    ),
    "seq_center": NextflowParameter(
        type=Optional[str],
        display_name="Sequencing Center",
        description="Sequencing center information for BAM read groups.",
    ),
    "save_unaligned": NextflowParameter(
        type=Optional[bool],
        display_name="Save Unaligned Reads",
        description="Save unaligned reads to the results directory.",
    ),
    "save_reference": NextflowParameter(
        type=Optional[bool],
        display_name="Save Reference Files",
        description="Save generated reference genome files and indices.",
        default=True,
    ),
    "genome_source": NextflowParameter(),
    "genome": NextflowParameter(
        type=Optional[str],
        display_name="iGenomes Reference",
        description="Name of iGenomes reference.",
    ),
    "fasta": NextflowParameter(
        type=Optional[LatchFile],
        display_name="Genome FASTA",
        description="Path to FASTA genome file.",
    ),
    "gtf": NextflowParameter(
        type=Optional[str],
        display_name="GTF Annotation",
        description="Path to reference GTF file.",
    ),
    "mature": NextflowParameter(
        type=Optional[str],
        display_name="Mature miRNA FASTA",
        description="FASTA file with mature miRNAs for interaction analysis.",
    ),
    "bowtie": NextflowParameter(
        type=Optional[str],
        display_name="Bowtie Index",
        description="Path to Bowtie index files.",
    ),
    "bowtie2": NextflowParameter(
        type=Optional[str],
        display_name="Bowtie2 Index",
        description="Path to Bowtie2 index files.",
    ),
    "bwa": NextflowParameter(
        type=Optional[str],
        display_name="BWA Index",
        description="Path to BWA index directory.",
    ),
    "hisat2": NextflowParameter(
        type=Optional[str],
        display_name="HISAT2 Index",
        description="Path to HISAT2 index directory.",
    ),
    "hisat2_build_memory": NextflowParameter(
        type=Optional[str],
        display_name="HISAT2 Build Memory",
        description="Memory required for HISAT2 index building.",
        default="200.GB",
    ),
    "segemehl": NextflowParameter(
        type=Optional[str],
        display_name="Segemehl Index",
        description="Path to Segemehl index file.",
    ),
    "star": NextflowParameter(
        type=Optional[str],
        display_name="STAR Index",
        description="Path to STAR index directory.",
    ),
    "skip_trimming": NextflowParameter(
        type=Optional[bool],
        display_name="Skip Trimming",
        description="Skip the adapter trimming step.",
    ),
    "save_trimmed": NextflowParameter(
        type=Optional[bool],
        display_name="Save Trimmed Reads",
        description="Save the trimmed FastQ files.",
    ),
    "skip_fastqc": NextflowParameter(
        type=Optional[bool],
        display_name="Skip FastQC",
        description="Skip FastQC quality control.",
    ),
    "clip_r1": NextflowParameter(
        type=Optional[int],
        display_name="Clip Read 1 (5')",
        description="Number of bases to remove from 5' end of read 1.",
    ),
    "clip_r2": NextflowParameter(
        type=Optional[int],
        display_name="Clip Read 2 (5')",
        description="Number of bases to remove from 5' end of read 2.",
    ),
    "three_prime_clip_r1": NextflowParameter(
        type=Optional[int],
        display_name="Clip Read 1 (3')",
        description="Number of bases to remove from 3' end of read 1 after trimming.",
    ),
    "three_prime_clip_r2": NextflowParameter(
        type=Optional[int],
        display_name="Clip Read 2 (3')",
        description="Number of bases to remove from 3' end of read 2 after trimming.",
    ),
    "trim_nextseq": NextflowParameter(
        type=Optional[int],
        display_name="NextSeq Trim Quality",
        description="Quality threshold for NextSeq trimming after poly-G removal.",
    ),
    "min_trimmed_reads": NextflowParameter(
        type=Optional[int],
        display_name="Min Trimmed Reads",
        description="Minimum number of reads required after trimming.",
        default=10000,
    ),
    "save_intermediates": NextflowParameter(
        type=Optional[bool],
        display_name="Save Intermediate Files",
        description="Save intermediate processing files.",
    ),
    "multiqc_methods_description": NextflowParameter(
        type=Optional[str],
        display_name="MultiQC Methods Description",
        description="Custom MultiQC yaml file with methods description HTML.",
    ),
}
