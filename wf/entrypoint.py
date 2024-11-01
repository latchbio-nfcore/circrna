import csv
import os
import shutil
import subprocess
import sys
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import List, Optional

import requests
import typing_extensions
from flytekit.core.annotation import FlyteAnnotation
from latch.executions import rename_current_execution, report_nextflow_used_storage
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata


@dataclass
class SampleSheet:
    sample: str
    fastq_1: LatchFile
    fastq_2: Optional[LatchFile]
    strandedness: Optional[str] = None


def custom_samplesheet_constructor(
    samples: List[SampleSheet], shared_dir: Path
) -> Path:
    samplesheet = shared_dir / "samplesheet.csv"

    columns = ["sample", "fastq_1", "fastq_2", "strandedness"]

    with open(samplesheet, "w") as f:
        writer = csv.DictWriter(f, columns, delimiter=",")
        writer.writeheader()

        for sample in samples:
            row_data = {
                "sample": sample.sample,
                "fastq_1": sample.fastq_1.remote_path,
                "fastq_2": sample.fastq_2.remote_path if sample.fastq_2 else "",
                "strandedness": sample.strandedness
                if sample.strandedness is not None
                else "auto",
            }
            writer.writerow(row_data)

    return samplesheet


@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize(run_name: str) -> str:
    rename_current_execution(run_name)
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage-ofs",
        headers=headers,
        json={
            "storage_expiration_hours": 0,
            "version": 2,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]


@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    run_name: str,
    pvc_name: str,
    input: LatchFile,
    outdir: LatchOutputDir,
    phenotype: Optional[LatchFile],
    annotation: Optional[LatchFile],
    email: Optional[str],
    multiqc_title: Optional[str],
    mirna_expression: Optional[LatchFile],
    seq_center: Optional[str],
    genome: Optional[str],
    fasta: Optional[LatchFile],
    gtf: Optional[str],
    mature: Optional[str],
    bowtie: Optional[str],
    bowtie2: Optional[str],
    bwa: Optional[str],
    hisat2: Optional[str],
    segemehl: Optional[str],
    star: Optional[str],
    clip_r1: Optional[int],
    clip_r2: Optional[int],
    three_prime_clip_r1: Optional[int],
    three_prime_clip_r2: Optional[int],
    trim_nextseq: Optional[int],
    multiqc_methods_description: Optional[str],
    tools: str,
    bsj_reads: Optional[int],
    max_shift: Optional[int],
    min_tools: Optional[int],
    min_samples: Optional[int],
    exon_boundary: Optional[int],
    quantification_tools: Optional[str],
    bootstrap_samples: Optional[int],
    mirna_min_sample_percentage: Optional[float],
    mirna_min_reads: Optional[int],
    mirna_correlation: Optional[str],
    mirna_tools: Optional[str],
    mirna_min_tools: Optional[int],
    sjdboverhang: Optional[int],
    chimJunctionOverhangMin: Optional[int],
    alignSJDBoverhangMin: Optional[int],
    limitSjdbInsertNsj: Optional[int],
    chimSegmentMin: Optional[int],
    seglen: Optional[int],
    min_intron: Optional[int],
    max_intron: Optional[int],
    min_map_len: Optional[int],
    min_fusion_distance: Optional[int],
    save_unaligned: bool,
    save_reference: bool,
    hisat2_build_memory: Optional[str],
    skip_trimming: bool,
    save_trimmed: bool,
    skip_fastqc: bool,
    min_trimmed_reads: Optional[int],
    save_intermediates: bool,
) -> None:
    shared_dir = Path("/nf-workdir")

    input_sheet = custom_samplesheet_constructor(input=input, shared_dir=shared_dir)

    exec_name = _get_execution_name()
    if exec_name is None:
        print("Failed to get execution name.")
        exec_name = "unknown"

    latch_log_dir = urljoins("latch:///your_log_dir/nf_nf_core_circrna", exec_name)
    print(f"Log directory: {latch_log_dir}")

    ignore_list = [
        "latch",
        ".latch",
        ".git",
        "nextflow",
        ".nextflow",
        "work",
        "results",
        "miniconda",
        "anaconda3",
        "mambaforge",
    ]

    shutil.copytree(
        Path("/root"),
        shared_dir,
        ignore=lambda src, names: ignore_list,
        ignore_dangling_symlinks=True,
        dirs_exist_ok=True,
    )

    cmd = [
        "/root/nextflow",
        "run",
        str(shared_dir / "main.nf"),
        "-work-dir",
        str(shared_dir),
        "-profile",
        "docker",
        "-c",
        "latch.config",
        "-resume",
        *get_flag("input", input_sheet),
        *get_flag("outdir", outdir),
        *get_flag("phenotype", phenotype),
        *get_flag("annotation", annotation),
        *get_flag("email", email),
        *get_flag("multiqc_title", multiqc_title),
        *get_flag("tools", tools),
        *get_flag("bsj_reads", bsj_reads),
        *get_flag("max_shift", max_shift),
        *get_flag("min_tools", min_tools),
        *get_flag("min_samples", min_samples),
        *get_flag("exon_boundary", exon_boundary),
        *get_flag("quantification_tools", quantification_tools),
        *get_flag("bootstrap_samples", bootstrap_samples),
        *get_flag("mirna_expression", mirna_expression),
        *get_flag("mirna_min_sample_percentage", mirna_min_sample_percentage),
        *get_flag("mirna_min_reads", mirna_min_reads),
        *get_flag("mirna_correlation", mirna_correlation),
        *get_flag("mirna_tools", mirna_tools),
        *get_flag("mirna_min_tools", mirna_min_tools),
        *get_flag("sjdboverhang", sjdboverhang),
        *get_flag("chimJunctionOverhangMin", chimJunctionOverhangMin),
        *get_flag("alignSJDBoverhangMin", alignSJDBoverhangMin),
        *get_flag("limitSjdbInsertNsj", limitSjdbInsertNsj),
        *get_flag("chimSegmentMin", chimSegmentMin),
        *get_flag("seglen", seglen),
        *get_flag("min_intron", min_intron),
        *get_flag("max_intron", max_intron),
        *get_flag("min_map_len", min_map_len),
        *get_flag("min_fusion_distance", min_fusion_distance),
        *get_flag("seq_center", seq_center),
        *get_flag("save_unaligned", save_unaligned),
        *get_flag("save_reference", save_reference),
        *get_flag("genome", genome),
        *get_flag("fasta", fasta),
        *get_flag("gtf", gtf),
        *get_flag("mature", mature),
        *get_flag("bowtie", bowtie),
        *get_flag("bowtie2", bowtie2),
        *get_flag("bwa", bwa),
        *get_flag("hisat2", hisat2),
        *get_flag("hisat2_build_memory", hisat2_build_memory),
        *get_flag("segemehl", segemehl),
        *get_flag("star", star),
        *get_flag("skip_trimming", skip_trimming),
        *get_flag("save_trimmed", save_trimmed),
        *get_flag("skip_fastqc", skip_fastqc),
        *get_flag("clip_r1", clip_r1),
        *get_flag("clip_r2", clip_r2),
        *get_flag("three_prime_clip_r1", three_prime_clip_r1),
        *get_flag("three_prime_clip_r2", three_prime_clip_r2),
        *get_flag("trim_nextseq", trim_nextseq),
        *get_flag("min_trimmed_reads", min_trimmed_reads),
        *get_flag("save_intermediates", save_intermediates),
        *get_flag("multiqc_methods_description", multiqc_methods_description),
    ]

    print("Launching Nextflow Runtime")
    print(" ".join(cmd))
    print(flush=True)

    failed = False
    try:
        env = {
            **os.environ,
            "NXF_ANSI_LOG": "false",
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms1536M -Xmx6144M -XX:ActiveProcessorCount=4",
            "NXF_DISABLE_CHECK_LATEST": "true",
            "NXF_ENABLE_VIRTUAL_THREADS": "false",
            "NXF_ENABLE_FS_SYNC": "true",
        }

        if False:
            env["LATCH_LOG_DIR"] = latch_log_dir

        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    except subprocess.CalledProcessError:
        failed = True
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            remote = LPath(urljoins(latch_log_dir, "nextflow.log"))
            print(f"Uploading .nextflow.log to {remote.path}")
            remote.upload_from(nextflow_log)

        print("Computing size of workdir... ", end="")
        try:
            result = subprocess.run(
                ["du", "-sb", str(shared_dir)],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=5 * 60,
            )

            size = int(result.stdout.split()[0])
            report_nextflow_used_storage(size)
            print(f"Done. Workdir size: {size / 1024 / 1024 / 1024: .2f} GiB")
        except subprocess.TimeoutExpired:
            print(
                "Failed to compute storage size: Operation timed out after 5 minutes."
            )
        except subprocess.CalledProcessError as e:
            print(f"Failed to compute storage size: {e.stderr}")
        except Exception as e:
            print(f"Failed to compute storage size: {e}")

    if failed:
        sys.exit(1)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_circrna(
    run_name: str,
    input: LatchFile,
    outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
    phenotype: Optional[LatchFile],
    annotation: Optional[LatchFile],
    email: Optional[str],
    multiqc_title: Optional[str],
    mirna_expression: Optional[LatchFile],
    seq_center: Optional[str],
    genome: Optional[str],
    fasta: Optional[LatchFile],
    gtf: Optional[str],
    mature: Optional[str],
    bowtie: Optional[str],
    bowtie2: Optional[str],
    bwa: Optional[str],
    hisat2: Optional[str],
    segemehl: Optional[str],
    star: Optional[str],
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
) -> None:
    """
    nf-core/circrna

    Sample Description
    """

    pvc_name: str = initialize(run_name=run_name)
    nextflow_runtime(
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
