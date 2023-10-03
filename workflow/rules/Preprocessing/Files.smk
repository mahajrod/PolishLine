localrules: create_paired_fastx_links, create_single_fastx_links, create_links_for_draft, transfer_input_files
ruleorder: extract_paired_fastq_from_bam > create_paired_fastx_links
ruleorder: extract_single_fastq_from_bam > create_single_fastx_links
ruleorder: create_links_for_draft > transfer_input_files
rule create_paired_fastx_links:
    priority: 1000
    input:
        forward_fastq=lambda wildcards: input_dir_path.resolve() / "{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                              wildcards.dataformat,
                                                                                              wildcards.pairprefix,
                                                                                              config["data"][wildcards.datatype]["forward_suffix"],
                                                                                              config["data"][wildcards.datatype]["input_extension"]),
        reverse_fastq=lambda wildcards: input_dir_path.resolve() / "{0}/{1}/{2}{3}{4}".format(wildcards.datatype,
                                                                                              wildcards.dataformat,
                                                                                              wildcards.pairprefix,
                                                                                              config["data"][wildcards.datatype]["reverse_suffix"],
                                                                                              config["data"][wildcards.datatype]["input_extension"]),
    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        forward_fastq=output_dict["data"] / ("{dataformat, fast[aq]}/{datatype, [^/]+}/raw/{pairprefix}_1%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("{dataformat, fast[aq]}/{datatype, [^/]+}/raw/{pairprefix}_2%s" % config["fastq_extension"])
    log:
        std=output_dict["log"] / "create_fastq_links.{datatype}.{dataformat}.{pairprefix}.log",
        cluster_log=output_dict["cluster_log"] / "create_fastq_links.{datatype}.{dataformat}.{pairprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_fastq_links.{datatype}.{dataformat}.{pairprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fastq_links.{datatype}.{dataformat}.{pairprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_fastq_links"),
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input.forward_fastq} {output.forward_fastq} > {log.std} 2>&1; "
         " ln -sf {input.reverse_fastq} {output.reverse_fastq} >> {log.std} 2>&1; "

rule create_single_fastx_links:
    priority: 1000
    input:
        fastq=lambda wildcards: input_dir_path.resolve() / "{0}/{1}/{2}{3}".format(wildcards.datatype,
                                                                                   wildcards.dataformat,
                                                                                   wildcards.fileprefix,
                                                                                   config["data"][wildcards.datatype]["input_extension"]),

    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        fastq=output_dict["data"] / ("{dataformat}/{datatype, [^/]+}/raw/{fileprefix}%s" % config["fastq_extension"]),
    log:
        std=output_dict["log"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_fastq_links"),
        cpus=parameters["threads"]["create_fastq_links"],
        time=parameters["time"]["create_fastq_links"],
        mem=parameters["memory_mb"]["create_fastq_links"],
    threads:
        parameters["threads"]["create_fastq_links"]
    shell:
         " ln -sf {input.fastq} {output.fastq} > {log.std} 2>&1; "

rule extract_single_fastq_from_bam:
    priority: 1000
    input:
        bam=lambda wildcards: input_dir_path.resolve() / "{0}/{1}/{2},bam".format(wildcards.datatype,
                                                                                   wildcards.dataformat,
                                                                                   wildcards.fileprefix),

    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        fastq=output_dict["data"] / ("{dataformat, mapped_bam|unmapped_bam}/{datatype, [^/]+}/raw/{fileprefix}%s" % config["fastq_extension"]),
    log:
        std=output_dict["log"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_fastq_from_bam"),
        cpus=parameters["threads"]["extract_fastq_from_bam"],
        time=parameters["time"]["extract_fastq_from_bam"],
        mem=parameters["memory_mb"]["extract_fastq_from_bam"],
    threads:
        parameters["threads"]["extract_fastq_from_bam"]
    shell:
         " samtools fastq -@ {threads} -s {output.fastq} {input.bam} > {log.std} 2>&1; "


rule extract_paired_fastq_from_bam:
    priority: 1000
    input:
        bam=lambda wildcards: input_dir_path.resolve() / "{0}/{1}/{2},bam".format(wildcards.datatype,
                                                                                  wildcards.dataformat,
                                                                                  wildcards.fileprefix),

    output:
        #directory(output_dict["data"] / "/fastq/{datatype}/raw"),
        forward_fastq=output_dict["data"] / ("{dataformat, mapped_bam|unmapped_bam}/{datatype, [^/]+}/raw/{fileprefix}_1%s" % config["fastq_extension"]),
        reverse_fastq=output_dict["data"] / ("{dataformat, mapped_bam|unmapped_bam}/{datatype, [^/]+}/raw/{fileprefix}_2%s" % config["fastq_extension"]),
    log:
        std=output_dict["log"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.log",
        cluster_log=output_dict["cluster_log"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "extract_single_fastq_from_bam.{datatype}.{dataformat}.{fileprefix}.cluster.err",
    benchmark:
        output_dict["benchmark"] / "create_fastq_links.{datatype}.{dataformat}.{fileprefix}.benchmark.txt",
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("extract_fastq_from_bam"),
        cpus=parameters["threads"]["extract_fastq_from_bam"],
        time=parameters["time"]["extract_fastq_from_bam"],
        mem=parameters["memory_mb"]["extract_fastq_from_bam"],
    threads:
        parameters["threads"]["extract_fastq_from_bam"]
    shell:
         " samtools fastq -@ {threads} -1 {output.forward_fastq} -2 {output.reverse_fastq} {input.bam} > {log.std} 2>&1; "

rule create_links_for_draft:
    input:
        lambda wildcards: input_dir_path.resolve() / "draft/fasta/{0}".format(draft_file_dict[wildcards.haplotype])
    output:
        out_dir_path / "draft_qc/{parameters}/{genome_prefix}.draft_qc.{haplotype}.fasta"
    log:
        ln=output_dict["log"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_draft.{genome_prefix}.{parameters}.draft_qc.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_links_for_draft"),
        cpus=parameters["threads"]["create_links_for_draft"],
        time=parameters["time"]["create_links_for_draft"],
        mem=parameters["memory_mb"]["create_links_for_draft"]
    threads: parameters["threads"]["create_links_for_draft"]

    shell:
        " ln -sf {input} {output} 2>{log.ln}; "

rule transfer_input_files: #
    input:
        fasta=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict[wildcards.stage]["prev_stage"],
                                                                                  wildcards.prev_stage_parameters,
                                                                                  wildcards.genome_prefix,
                                                                                  wildcards.haplotype),
        #len=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.len" % (stage_dict["curation"]["prev_stage"],
        #                                                                                           stage_dict["curation"]["prev_stage"])),
        #fai=out_dir_path / ("%s/{prev_stage_parameters}/{genome_prefix}.%s.{haplotype}.fasta.fai" % (stage_dict["curation"]["prev_stage"],
        #                                                                                               stage_dict["curation"]["prev_stage"])),
        #bed=get_hic_bed_file if not config["skip_higlass"] else []
    output:
        fasta=out_dir_path / "{stage}/{prev_stage_parameters}..{current_stage_parameters}/{genome_prefix}.{stage}.{haplotype, [^.]+}.fasta",
        #len=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.len",
        #fai=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/scaffolds/{genome_prefix}.input.{haplotype}.fasta.fai",
        #bed=out_dir_path / "curation/{prev_stage_parameters}..{curation_parameters}/{haplotype, [^.]+}/input/{genome_prefix}.input.{haplotype}.hic.bed" if not config["skip_higlass"] else [],
    log:
        cp=output_dict["log"]  / "create_input_files.{stage}.{prev_stage_parameters}..{current_stage_parameters}.{genome_prefix}.{haplotype}.cp.log",
        cluster_log=output_dict["cluster_log"] / "create_input_files.{stage}.{prev_stage_parameters}..{current_stage_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_input_files.{stage}.{prev_stage_parameters}..{current_stage_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "create_phasing_files.{stage}.{prev_stage_parameters}..{current_stage_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("create_input_files"),
        cpus=parameters["threads"]["create_input_files"],
        time=parameters["time"]["create_input_files"],
        mem=parameters["memory_mb"]["create_input_files"]
    threads: parameters["threads"]["create_input_files"]

    shell:
        " cp -f `realpath -s {input.fasta}` {output.fasta} > {log.cp} 2>&1; "
        #" cp -f `realpath -s {input.fai}` {output.fai} >> {log.cp} 2>&1; "
        #" cp -f `realpath -s {input.len}` {output.len} >> {log.cp} 2>&1; "
