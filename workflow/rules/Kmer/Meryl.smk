#ruleorder: meryl_pe > create_fastq_links
rule meryl:
    input:
        lambda wildcards: output_dict["data"] / "{0}/{1}/{2}/{3}{4}".format(config["data"][wildcards.datatype]["converted_format"],
                                                                            wildcards.datatype,
                                                                            wildcards.stage,
                                                                            wildcards.fileprefix,
                                                                            config["data"][wildcards.datatype]["converted_extension"])
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length, [^./]+}.meryl.{fileprefix, (?!^histo$)}") #, (?!^histo$)
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{fileprefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("meryl"),
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=partial(set_mem_limit,
                    default_mem=parameters["memory_mb"]["meryl"],
                    multiplicator=1),
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"


rule meryl_pe:
    input:
        forward_fastq=lambda wildcards: output_dict["data"] / "{0}/{1}/{2}/{3}{4}{5}".format(config["data"][wildcards.datatype]["converted_format"],
                                                                                             wildcards.datatype,
                                                                                             wildcards.stage,
                                                                                             wildcards.pairprefix,
                                                                                             config["data"][wildcards.datatype]["converted_forward_suffix"],
                                                                                             config["data"][wildcards.datatype]["converted_extension"]),
        reverse_fastq=lambda wildcards: output_dict["data"] / "{0}/{1}/{2}/{3}{4}{5}".format(config["data"][wildcards.datatype]["converted_format"],
                                                                                             wildcards.datatype,
                                                                                             wildcards.stage,
                                                                                             wildcards.pairprefix,
                                                                                             config["data"][wildcards.datatype]["converted_reverse_suffix"],
                                                                                             config["data"][wildcards.datatype]["converted_extension"]),
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.{pairprefix}") # , (?!^histo$)
    log:
        std=output_dict["log"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.log",
        cluster_log=output_dict["cluster_log"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "meryl.{datatype}.{stage}.{pairprefix}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("meryl_pe"),
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=partial(set_mem_limit,
                    default_mem=parameters["memory_mb"]["meryl"],
                    multiplicator=1),
        kmer_counter=1
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl k={wildcards.kmer_length} threads={threads} memory={resources.mem}m count "
         " output {output.db_dir} {input} 1>{log.std} 2>&1;"

rule merge_meryl:
    input:
        lambda wildcards:
            expand(output_dict["kmer"] / ("%s/%s/%s.%s.%s.meryl.{fileprefix}" % (wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 wildcards.datatype,
                                                                                 wildcards.stage,
                                                                                 wildcards.kmer_length,)),
                   #fileprefix=input_file_prefix_dict[wildcards.datatype] if datatype_format_dict[wildcards.datatype] == "fastq" else input_fasta_file_prefix_dict[wildcards.datatype],
                   fileprefix=config["data"][wildcards.datatype]["converted_fileprefix_list"],
                   allow_missing=True,)  if wildcards.datatype not in config["paired_fastq_based_data"] else \
            expand(rules.meryl_pe.output,
                   pairprefix=config["data"][wildcards.datatype]["pairprefix_list"], # input_pairprefix_dict[wildcards.datatype],
                   allow_missing=True,)
    output:
        db_dir=directory(output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl"),
        histo=output_dict["kmer"] / "{datatype}/{stage}/{datatype}.{stage}.{kmer_length}.meryl.histo"

    log:
        count_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.count.log",
        histo_log=output_dict["log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.histo.log",
        cluster_log=output_dict["cluster_log"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.cluster.err"
    benchmark:
        output_dict["benchmark"] / "merge_meryl.{datatype}.{stage}.{kmer_length}.benchmark.txt"
    conda:
        config["conda"]["common"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["common"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("merge_meryl"),
        cpus=parameters["threads"]["meryl"],
        time=parameters["time"]["meryl"],
        mem=partial(set_mem_limit,
                    default_mem=parameters["memory_mb"]["meryl"],
                    multiplicator=1),
    threads:
        parameters["threads"]["meryl"]
    shell:
         " meryl threads={threads} memory={resources.mem}m"
         " union-sum output {output.db_dir} {input} 1>{log.count_log} 2>&1;"
         " meryl threads={threads} memory={resources.mem}m "
         " histogram {output.db_dir} > {output.histo} 2>{log.histo_log}"
