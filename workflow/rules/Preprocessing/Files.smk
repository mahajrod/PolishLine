localrules: create_paired_fastx_links, create_single_fastx_links, create_links_for_draft,
ruleorder: extract_paired_fastq_from_bam > create_paired_fastx_links
ruleorder: extract_single_fastq_from_bam > create_single_fastx_links

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

"""
rule create_links_for_reference:
    input:
        fasta=lambda wildcards: input_reference_filedict[wildcards.ref_name]["fasta"].resolve(),
        syn=lambda wildcards: input_reference_filedict[wildcards.ref_name]["syn"].resolve(),
        whitelist=lambda wildcards: input_reference_filedict[wildcards.ref_name]["whitelist"].resolve(),
        orderlist=lambda wildcards: input_reference_filedict[wildcards.ref_name]["orderlist"].resolve(),
    output:
        fasta=out_dir_path / "data/reference/{ref_name}/{ref_name}.softmasked.fasta",
        syn=out_dir_path / "data/reference/{ref_name}/{ref_name}.syn",
        whitelist=out_dir_path / "data/reference/{ref_name}/{ref_name}.whitelist",
        orderlist=out_dir_path / "data/reference/{ref_name}/{ref_name}.orderlist",
    log:
        ln=output_dict["log"]  / "create_links_for_reference.{ref_name}.ln.log",
        cluster_log=output_dict["cluster_log"] / "create_links_for_reference.{ref_name}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "create_links_for_reference.{ref_name}.err"
    benchmark:
        output_dict["benchmark"]  / "create_links_for_reference.{ref_name}.benchmark.txt"
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
        " ln -sf {input.fasta} {output.fasta} 2>{log.ln}; "
        " ln -sf {input.syn} {output.syn} 2>>{log.ln}; "
        " ln -sf {input.whitelist} {output.whitelist} 2>>{log.ln}; "
        " ln -sf {input.orderlist} {output.orderlist} 2>>{log.ln}; "
"""