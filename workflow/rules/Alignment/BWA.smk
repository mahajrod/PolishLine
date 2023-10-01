

rule bwa_map_pe: #
    input:
        index=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{genome_prefix}.{stage}.{haplotype}.fasta.ann",
        reference=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{genome_prefix}.{stage}.{haplotype}.fasta",
        forward_fastq=lambda wildcards: out_dir_path / "{0}/{1}/fastq/{2}/{3}/{4}/{5}_1.fastq.gz".format(stage_dict[wildcards.stage]["prev_stage"],
                                                                                                              wildcards.prev_parameters,
                                                                                                              wildcards.haplotype,
                                                                                                              wildcards.phasing_kmer_length,
                                                                                                              wildcards.datatype,
                                                                                                              wildcards.pairprefix),# TODO: change to forward_suffix
        reverse_fastq=lambda wildcards: out_dir_path / "{0}/{1}/fastq/{2}/{3}/{4}/{5}_2.fastq.gz".format(stage_dict[wildcards.stage]["prev_stage"],
                                                                                                              wildcards.prev_parameters,
                                                                                                              wildcards.haplotype,
                                                                                                              wildcards.phasing_kmer_length,
                                                                                                              wildcards.datatype,
                                                                                                              wildcards.pairprefix),# TODO: change to forward_suffix
    output:
        #bam=out_dir_path  / "{assembly_stage}/{parameters}/{haplotype}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.{fileprefix}.bwa.bam"
        bam=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{haplotype, [^.]+}/{aligner, bwa-mem2|bwa}/{datatype}/{phasing_kmer_length, [^./]+}/{genome_prefix}.{stage}.{datatype}.{phasing_kmer_length}.{haplotype}.{pairprefix, [^/]+}.{aligner}.bam"
    params:
        id=config["genome_prefix"],
    log:
        map=output_dict["log"]  / "bwa_map.{aligner}.{stage}.{prev_parameters}..{alignment_parameters}.{datatype}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.map.log",
        sort=output_dict["log"]  / "bwa_map.{aligner}.{stage}.{prev_parameters}..{alignment_parameters}.{datatype}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.sort.log",
        #filter=output_dict["log"]  / "bwa_map.{assembly_stage}.{parameters}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.filter.log",
        cluster_log=output_dict["cluster_log"] / "bwa_map.{aligner}.{stage}.{prev_parameters}..{alignment_parameters}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bwa_map.{aligner}.{stage}.{prev_parameters}..{alignment_parameters}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.{pairprefix}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bwa_map.{aligner}.{stage}.{prev_parameters}..{alignment_parameters}.{datatype}.{genome_prefix}.{haplotype}.{phasing_kmer_length}.{pairprefix}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bwa_map"),
        cpus=parameters["threads"]["bwa_map"] ,
        time=parameters["time"]["bwa_map"],
        mem=parameters["memory_mb"]["bwa_map"]
    threads: parameters["threads"]["bwa_map"]
    shell:
        " {wildcards.aligner} mem -SP5M -t {threads} -R  \'@RG\\tID:{params.id}\\tPU:x\\tSM:{params.id}\\tPL:illumina\\tLB:x\' "
        " {input.reference} {input.forward_fastq} {input.reverse_fastq} 2>{log.map} | samtools view -Sb - > {output.bam} 2>{log.sort} "

rule bam_merge_pe_files:
    input:
        bams=lambda wildcards: expand(out_dir_path / "{0}/{1}..{2}/{3}/{4}/{5}/{6}/{7}.{0}.{5}.{6}.{3}.{{pairprefix}}.{4}.bam".format(wildcards.stage,
                                                                                                                                      wildcards.prev_parameters,
                                                                                                                                      wildcards.alignment_parameters,
                                                                                                                                      wildcards.haplotype,
                                                                                                                                      wildcards.aligner,
                                                                                                                                      wildcards.datatype,
                                                                                                                                      wildcards.phasing_kmer_length,
                                                                                                                                      wildcards.genome_prefix),
                    allow_missing=True,
                    pairprefix=config["data"][wildcards.datatype]["pairprefix_list"]), #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        reference_fai=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{genome_prefix}.{stage}.{haplotype}.fasta.fai",
        reference=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{genome_prefix}.{stage}.{haplotype}.fasta"
    output:
        bam=out_dir_path / "{stage, [^/]+}/{prev_parameters}..{alignment_parameters}/{haplotype, [^./]+}/{aligner, [^/]+}/{datatype, [^/]+}/{phasing_kmer_length, [^./]+}/{genome_prefix}.{stage}.{datatype}.{phasing_kmer_length}.{haplotype}.{aligner}.bam"
    params:
        sort_threads=parameters["threads"]["samtools_sort"]
    log:
        std=output_dict["log"] / "bam_merge_files.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.log",
        cluster_log=output_dict["cluster_log"] / "bam_merge_files.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "bam_merge_files.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "bam_merge_files.{stage}.{prev_parameters}..{alignment_parameters}.{genome_prefix}.{aligner}.{datatype}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("bwa_merge_files"),
        cpus=parameters["threads"]["samtools_sort"] ,
        time=parameters["time"]["samtools_sort"],
        mem=parameters["memory_mb"]["samtools_sort"]
    threads: parameters["threads"]["samtools_sort"]
    shell:
        " samtools merge -@ {params.sort_threads} -o {output.bam} {input.bams} 1>{log.std} 2>&1"

rule markdup:
    input:
        bam=out_dir_path / "{stage}/{prev_parameters}..{alignment_parameters}/{haplotype}/{aligner}/{datatype}/{phasing_kmer_length}/{genome_prefix}.{stage}.{datatype}.{phasing_kmer_length}.{haplotype}.{aligner}.bam"
    output:
        bam=out_dir_path / "{stage, [^/]+}/{prev_parameters}..{alignment_parameters}/{haplotype, [^./]+}/{aligner, [^/]+}/{datatype, [^/]+}/{phasing_kmer_length, [^./]+}/{genome_prefix}.{stage}.{datatype}.{phasing_kmer_length}.{haplotype}.{aligner}.markdup.bam"
        #bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
    params:
        sort_threads=parameters["threads"]["samtools_sort"],
        collate_threads=parameters["threads"]["samtools_collate"],
        fixmate_threads=parameters["threads"]["samtools_fixmate"],
        markdup_threads=parameters["threads"]["samtools_markdup"],
        sort_per_thread=parameters["memory_mb"]["samtools_sort"]
    log:
        collate=output_dict["log"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.collate.log",
        fixmate=output_dict["log"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.fixmate.log",
        sort=output_dict["log"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.sort.log",
        markdup=output_dict["log"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.markdup.log",
        cluster_log=output_dict["cluster_log"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "markdup.{stage}.{prev_parameters}..{alignment_parameters}.{aligner}.{datatype}.{genome_prefix}.{phasing_kmer_length}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("markdup"),
        cpus=parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"],
        time=parameters["time"]["markdup"],
        mem=partial(set_mem_limit,
                    default_mem=10000 + parameters["memory_mb"]["samtools_collate"] + parameters["memory_mb"]["samtools_fixmate"] + parameters["memory_mb"]["samtools_markdup"] + parameters["memory_mb"]["samtools_sort"] * parameters["threads"]["samtools_sort"],
                    multiplicator=1)
    threads: parameters["threads"]["samtools_sort"] + parameters["threads"]["samtools_collate"] + parameters["threads"]["samtools_fixmate"] + parameters["threads"]["samtools_markdup"]
    shell:
        " TMP_DIR=`dirname {output.bam}`; "
        " samtools collate -T ${{TMP_DIR}}/tmp.collate  -@ {params.collate_threads}  -O {input.bam} 2>{log.collate} | "
        " samtools fixmate -@ {params.fixmate_threads} -m - -  2>{log.fixmate} | "
        " samtools sort -T ${{TMP_DIR}}/tmp.sort -@ {params.sort_threads} -m {params.sort_per_thread}M 2>{log.sort} | "
        " samtools markdup -@ {params.markdup_threads} - {output.bam} > {log.markdup} 2>&1; "
