ruleorder: pilon > transfer_input_files
rule pilon:
    input:
        bam=lambda wildcards: out_dir_path / "{0}/{1}/{2}/{3}/{4}/{5}/{6}.{0}.{4}.{5}.{2}.{3}.markdup.bam".format(stage_dict[wildcards.stage]["prev_stage"],
                                                                                                                   wildcards.prev_parameters,
                                                                                                                   wildcards.haplotype,
                                                                                                                   stage_dict[stage_dict[wildcards.stage]["prev_stage"]]["parameters"][wildcards.prev_parameters]["coretool"],#wildcards.aligner,
                                                                                                                   stage_dict[stage_dict[wildcards.stage]["prev_stage"]]["parameters"][wildcards.prev_parameters]["option_set"]["datatype"],
                                                                                                                   stage_dict[stage_dict[wildcards.stage]["prev_stage"]]["parameters"][wildcards.prev_parameters]["option_set"]["phasing_kmer_length"],# wildcards.phasing_kmer_length,
                                                                                                                   wildcards.genome_prefix),
        reference=lambda wildcards: out_dir_path / "{0}/{1}/{2}.{0}.{3}.fasta".format(stage_dict[wildcards.stage]["prev_stage"],
                                                                                      wildcards.prev_parameters,
                                                                                      wildcards.genome_prefix,
                                                                                      wildcards.haplotype)

    output:
        fasta_alias=out_dir_path / "{stage, [^/]+}/{prev_parameters}..{polish_parameters}/{genome_prefix}.{stage}.{haplotype, [^./]+}.fasta",
        fasta=out_dir_path / "{stage, [^/]+}/{prev_parameters}..{polish_parameters}/{haplotype, [^./]+}/{genome_prefix}.{stage}.{haplotype}.fasta"

        #bai=out_dir_path / "{assembly_stage}/{parameters}/{haplotype, [^.]+}/alignment/{phasing_kmer_length}/{genome_prefix}.{assembly_stage}.{phasing_kmer_length}.{haplotype}.rmdup.bam.bai",
    params:
        fix_list= lambda wildcards: " --fix {0}".format(",".join(stage_dict[wildcards.stage]["parameters"][wildcards.prev_parameters + ".." + wildcards.polish_parameters]["option_set"]["fix_list"])),
        data_ploidy=lambda wildcards: "" if stage_dict[stage_dict[wildcards.stage]["prev_stage"]]["parameters"][wildcards.prev_parameters]["option_set"]["phasing_kmer_length"] != "NA" else  \
                "" if stage_dict[wildcards.stage]["parameters"][wildcards.prev_parameters + ".." + wildcards.polish_parameters]["option_set"]["assembly_ploidy"] == 1 else "--diploid",
    log:
        std=output_dict["log"] / "pilon.{stage}.{prev_parameters}..{polish_parameters}.{genome_prefix}.{haplotype}.std.log",
        ln=output_dict["log"] / "pilon.{stage}.{prev_parameters}..{polish_parameters}.{genome_prefix}.{haplotype}.ln.log",
        cluster_log=output_dict["cluster_log"] / "pilon.{stage}.{prev_parameters}..{polish_parameters}.{genome_prefix}.{haplotype}.cluster.log",
        cluster_err=output_dict["cluster_error"] / "pilon.{stage}.{prev_parameters}..{polish_parameters}.{genome_prefix}.{haplotype}.cluster.err"
    benchmark:
        output_dict["benchmark"]  / "pilon.{stage}.{prev_parameters}..{polish_parameters}.{genome_prefix}.{haplotype}.benchmark.txt"
    conda:
        config["conda"]["pilon"]["name"] if config["use_existing_envs"] else ("../../../%s" % config["conda"]["pilon"]["yaml"])
    resources:
        queue=config["queue"]["cpu"],
        node_options=parse_node_list("pilon"),
        cpus=parameters["threads"]["pilon"],
        time=parameters["time"]["pilon"],
        mem=partial(set_mem_limit,
                    default_mem=parameters["memory_mb"]["pilon"],
                    multiplicator=0.8)
    threads: parameters["threads"]["pilon"]
    shell:
        " OUTPUT_DIR=`dirname {output.fasta}`; "
        " OUTPUT_PREFIX=`basename {output.fasta}`; "
        " OUTPUT_PREFIX=${{OUTPUT_PREFIX%.fasta}}; "
        " pilon -Xmx{resources.mem}m --changes --vcf --tracks "
        " {params.data_ploidy} {params.fix_list} --genome {input.reference} --frags {input.bam} "
        " --output ${{OUTPUT_PREFIX}} --outdir ${{OUTPUT_DIR}} > {log.std} 2>&1; "
        " ln -sf {wildcards.haplotype}/`basename {output.fasta}` {output.fasta_alias} > {log.ln} 2>&1; "
