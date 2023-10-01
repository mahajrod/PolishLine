import os
import sys

import yaml
#import logging
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath

import pandas as pd

#---- Read config files ----
#-------- Read core config file --------
with open(config["main_config_file"], "r") as core_yaml_fd:
    config.update(yaml.safe_load(core_yaml_fd))
#---------------------------------------
#-------- Read resources config files --------
for resource, res_datatype in zip(["threads", "memory_mb", "time"], [int, int, str]):
    resource_df = pd.read_csv(config["resources"][resource], sep="\t", header=0, index_col=0)
    for config_label in resource_df.columns:
        config["parameters"][config_label][resource] = resource_df[config_label].to_dict(OrderedDict)

#---------------------------------------------

#---------------------------

#---- Include sections for functions ----
include: "workflow/functions/option_parsing.py"
include: "workflow/functions/general_parsing.py"
#----------------------------------------


#-- Initialization of path variables from config file --
#logging.info("Initialization of path variables...")
#---- Initialization of path variables for input----
input_dir_path = Path(config["input_dir"])

#----
#---- Initialization of path variables for output ----
out_dir_path = Path(config["out_dir"])
output_dict = {}

for first_level_sub_dir in config["first_level_subdir_list"]:
    output_dict[first_level_sub_dir] = out_dir_path / first_level_sub_dir

#----
#---- Initialization path variables for resources ----
#----

#---- Checking input files ----
#logging.info("Checking input files...")

data_types = config["data_types"].split(",")
for d_type in data_types:
    if d_type not in config["allowed_data_types"]:
        #logging.error("Unknown data type: {0}".format(d_type))
        raise ValueError("ERROR!!! Unknown data type: {0}".format(d_type))

config["data"] = {datatype: {} for datatype in data_types}
fastq_based_data_type_set = set()
fasta_based_data_type_set = set()
for datatype in data_types:
    datatype_dir = input_dir_path / datatype
    input = detect_input_type(datatype, input_dir_path / datatype)
    input_format = list(input.keys())[0]

    input_extension = list(input[input_format].keys())[0]
    config["data"][datatype] = {"input_dir": input_dir_path / datatype,
                                "input_format": input_format,
                                "input_extension": input_extension,
                                "file_list": deepcopy(input[input_format][input_extension]),
                                "fileprefix_list": list(map(lambda s: str(s.name)[:-len(input_extension)],
                                                            input[input_format][input_extension])),
                                "forward_suffix": None,
                                "reverse_suffix": None,
                                "pairprefix_list": [],
                                "converted_fileprefix_list": [],
                                "converted_forward_suffix": config["allowed_data"][datatype][input_format]["converted_forward_suffix"],
                                "converted_reverse_suffix": config["allowed_data"][datatype][input_format]["converted_reverse_suffix"],
                                "converted_format": config["allowed_data"][datatype][input_format]["converted_format"],
                                "converted_extension": config["allowed_data"][datatype][input_format]["converted_extension"]
                                }

    if config["data"][datatype]["converted_format"] == "fasta":
        fasta_based_data_type_set.add(datatype)
    elif config["data"][datatype]["converted_format"] == "fastq":
        fastq_based_data_type_set.add(datatype)

    if config["allowed_data"][datatype][input_format]["paired"]:
        # check filenames of paired data
        if (len(config["data"][datatype]["file_list"]) % 2) != 0:
            raise ValueError("ERROR!!! {0} fastq files seems to be unpaired or misrecognized".format(datatype))
        for forward, reverse in zip(config["data"][datatype]["file_list"][::2],
                                    config["data"][datatype]["file_list"][1::2]):
            if p_distance(str(forward), str(reverse), len(str(forward))) > 1:
                raise ValueError("ERROR!!! Forward and reverse read files differs by more than one symbol:\n\t{0}\n\t{1}".format(str(forward),
                                                                                                                                 str(reverse)))
        config["data"][datatype]["forward_suffix"] = set()
        config["data"][datatype]["reverse_suffix"] = set()
        # detect pairprefix, forward_and_reverse_suffixes for paired data
        for forward_prefix, reverse_prefix in zip(config["data"][datatype]["fileprefix_list"][::2],
                                                  config["data"][datatype]["fileprefix_list"][1::2]):
            common_prefix, forward_suffix, reverse_suffix = get_common_prefix_ans_suffixes(forward_prefix, reverse_prefix)
            config["data"][datatype]["pairprefix_list"].append(common_prefix)
            config["data"][datatype]["forward_suffix"].add(forward_suffix)
            config["data"][datatype]["reverse_suffix"].add(reverse_suffix)
        if (len(config["data"][datatype]["forward_suffix"]) > 1) or (len(config["data"][datatype]["reverse_suffix"]) > 1):
            raise ValueError("ERROR!!! Multiple different suffixes in filenames of %s data!" % d_type)

        config["data"][datatype]["forward_suffix"] = list(config["data"][datatype]["forward_suffix"])[0]
        config["data"][datatype]["reverse_suffix"] = list(config["data"][datatype]["reverse_suffix"])[0]

        for pairprefix in config["data"][datatype]["pairprefix_list"]:
            config["data"][datatype]["converted_fileprefix_list"].append(pairprefix + config["data"][datatype]["forward_suffix"])
            config["data"][datatype]["converted_fileprefix_list"].append(pairprefix + config["data"][datatype]["reverse_suffix"])
    else:
        config["data"][datatype]["converted_fileprefix_list"] = deepcopy(config["data"][datatype]["fileprefix_list"])

fastqc_data_type_set = fastq_based_data_type_set & set(config["fastqc_data_types"])
long_read_data_type_set = set(data_types) & set(config["long_read_data"])
genome_size_estimation_data_type_set = set(config["genome_size_estimation_data"]) & fastq_based_data_type_set & set(data_types)
coverage_track_data_type_set = set(data_types) & set(config["coverage_track_data"])
#---- Initialize tool parameters ----
#logging.info("Initializing tool parameters...")

if config["parameter_set"] not in config["parameters"]:
    raise ValueError("Error!!! Unknown set of tool parameters: {0}".format(config["parameter_set"]))

copy_absent_entries(config["parameters"]["default"], config["parameters"][config["parameter_set"]]) # set default values for options absent in  "parameter_set"

for key in list(config["parameters"].keys()): # remove unused sets of parameters
    if key != config["parameter_set"]:
        config["parameters"].pop(key)

parameters = config["parameters"][config["parameter_set"]] # short alias for used set of parameters
#print(parameters)

for tool in config["other_tool_option_sets"]: # select active set of option for tools other than coretools
    parameters["tool_options"][tool] = parameters["tool_options"][tool][config["other_tool_option_sets"][tool]]

#Kraken scan datatype
kraken_scan_data_type_set = set(data_types) & set(config["kraken_scan_data"])

#----
#---- Configure stages ----
config["stage_list"] = []

# Select configuration and combine stages from all mega_stages in a single list without nesting
if config["mode"] == "preprocessing":
    mega_stage_list = ["preprocessing"]
elif config["mode"] == "qc":
    mega_stage_list = ["preprocessing", "qc"]
elif config["mode"] == "phase":
    mega_stage_list = ["preprocessing", "qc", "phase"]
elif config["mode"] == "polish":
    mega_stage_list = ["preprocessing", "qc", "polish"]
else:
    raise ValueError("ERROR!!! Unknown mode: %s" % config["mode"])

for mega_stage in mega_stage_list:
    custom_megastage_entry = "custom_" + mega_stage + "_stages"
    if (custom_megastage_entry in config) and (config[custom_megastage_entry]):
        config["stage_list"].append(config[custom_megastage_entry])
    else:
        config["stage_list"] += config["allowed_stage_list"][mega_stage][config[mega_stage + "_mode"]]

stage_dict = OrderedDict()
for stage, stage_index in zip(config["stage_list"], range(0, len(config["stage_list"]))):
    stage_dict[stage] = OrderedDict()
    stage_dict[stage]["prev_stage"] = None if stage_index == 0 else config["stage_list"][stage_index-1]
print(stage_dict)
#----
#---- Save configuration and input files ----
final_config_yaml = output_dict["config"] / "config.final.yaml"
final_input_yaml = output_dict["config"] / "input.final.yaml"

os.makedirs(output_dict["config"], exist_ok=True)

with open(final_config_yaml, 'w') as final_config_fd, open(final_input_yaml, 'w') as final_input_fd:
    yaml.dump(convert_posixpath2str_in_dict(config), final_config_fd, default_flow_style=False, sort_keys=False)
    #yaml.dump(convert_posixpath2str_in_dict(input_dict), final_input_fd, default_flow_style=False, sort_keys=False)

#-------------------------------------------
localrules: all
#ruleorder: create_fastq_links > fastqc

results_dict = {}

haplotype_list = ["hap{0}".format(i) for i in range(1, config["ploidy"] + 1)] # TODO: obsolete: remove and fix issues
primary_haplotype = "hap1" # TODO: obsolete: remove and fix issues

results_list = []

for conda_env in config["conda"]:
    if "pip" in config["conda"][conda_env]:
        results_list += [expand("results/config/pip.{conda_env}.requirements",
                                conda_env=[conda_env])]

#---- Create output filelist ----
if "check_reads" in config["stage_list"]:
    results_list += [
                     final_config_yaml,
                     final_input_yaml
                     ]

if "check_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement


if ("read_qc" in config["stage_list"]) and ("read_qc" not in  config["skip_stage_dict"]):
    results_list += [*[expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=config["data"][datatype]["converted_fileprefix_list"],) for dat_type in fastqc_data_type_set ],
                      expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                             datatype=fastqc_data_type_set ,
                             stage=["raw",]),
                      *[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=config["data"][datatype]["converted_fileprefix_list"],) for dat_type in long_read_data_type_set],
                    *[expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                               datatype=[dat_type, ],
                               stage=["raw", ],
                               fileprefix=config["data"][datatype]["converted_fileprefix_list"],) for dat_type in long_read_data_type_set],
                     ]

if "draft_qc" in config["stage_list"]:
    draft_file_dict = get_input_assemblies(input_dir_path / "draft/fasta", config["ploidy"], config["assembly_fasta_extension"])
    stage_dict["draft_qc"]["parameters"] = {}

    for qcer in config["stage_coretools"]["draft_qc"]["default"]:
        for option_set in config["coretool_option_sets"][qcer]:
            parameters_label="{0}_{1}".format(qcer, option_set)
            stage_dict["draft_qc"]["parameters"][parameters_label] = {}
            stage_dict["draft_qc"]["parameters"][parameters_label]["qcer"] = qcer
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"] = {}
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["draft_qc"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
            stage_dict["draft_qc"]["parameters"][parameters_label]["option_set_group"] = None

    parameters_list = list(stage_dict["draft_qc"]["parameters"].keys())

    results_list += [expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["draft_qc"],),
                     *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                                assembly_stage=["draft_qc"],
                                parameters=[parameters_label],
                                genome_prefix=[config["genome_prefix"], ],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"]
                                ) for parameters_label in parameters_list],
                     ]
    if "busco" not in config["skip_stage_dict"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["draft_qc"],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["draft_qc"],
                                haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         ]
    #TODO: remove after debugging
    """
    results_list += [ *[expand(out_dir_path / "{stage}/{parameters}/kmer/{genome_prefix}.{stage}.{haplotype}.{assembly_kmer_length}",
                               stage=["draft_qc"],
                              parameters=[parameters_label],
                              genome_prefix=[config["genome_prefix"], ],
                              haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                              assembly_kmer_length=[31]) for parameters_label in parameters_list],
                      *[expand(out_dir_path / "{stage}/{parameters}/fasta/{haplotype}/{assembly_kmer_length}/{datatype}/{fileprefix}.fasta.gz",
                                stage=["draft_qc"],
                              parameters=[parameters_label],
                                datatype=[config["gap_closing_datatype"]],
                                fileprefix=input_file_prefix_dict[config["gap_closing_datatype"]] if datatype_format_dict[config["gap_closing_datatype"]] == "fastq" else input_fasta_file_prefix_dict[config["gap_closing_datatype"]],
                              genome_prefix=[config["genome_prefix"], ],
                              haplotype=stage_dict["draft_qc"]["parameters"][parameters_label]["haplotype_list"],
                              assembly_kmer_length=[31]) for parameters_label in parameters_list],

                      ]
    """


if ("filter_reads" in config["stage_list"]) and ("filter_reads" not in  config["skip_stage_dict"]):
    results_list += [expand(output_dict["data"] / ("fastq/hifi/filtered/{fileprefix}%s" % config["fastq_extension"]),
                            fileprefix=input_file_prefix_dict["hifi"]) if "hifi" in fastq_based_data_type_set else [],
                    expand(output_dict["qc"] / "fastqc/{datatype}/{stage}/{fileprefix}_fastqc.zip",
                           datatype=["hifi", ],
                           stage=["filtered", ],
                           fileprefix=input_file_prefix_dict["hifi"],
                           ) if "hifi" in fastq_based_data_type_set else [],
                    expand(output_dict["qc"] / "multiqc/{datatype}/{stage}/multiqc.{datatype}.{stage}.report.html",
                           datatype=["hifi"],
                           stage=["filtered",]) if "hifi" in fastq_based_data_type_set else [],
                    *[[expand(output_dict["kmer"] / "{datatype}/{stage}/genomescope/{genome_prefix}.{datatype}.{stage}.{kmer_length}.{kmer_tool}.genomescope.parameters",
                           datatype=[dat_type,],
                           genome_prefix=[config["genome_prefix"], ],
                           stage=["filtered",],
                           kmer_tool=[kmer_tool,],
                           kmer_length=parameters["tool_options"][kmer_tool][dat_type]["kmer_length"],
                           ) for kmer_tool in config["kmer_counter_list"] ]  for dat_type in genome_size_estimation_data_type_set],
                    ]
    if "nanoqc" not in config["skip_stage_dict"]:
        results_list += [

                        *[expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                                   datatype=[dat_type, ],
                                   stage=["filtered", ],
                                   fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                        expand(output_dict["qc"] / "nanoqc/{datatype}/{stage}/{fileprefix}",
                                   datatype=["nanopore", ],
                                   stage=["trimmed", ],
                                   fileprefix=input_file_prefix_dict["nanopore"],) if "nanopore" in long_read_data_type_set else [],
                        ]
    if "nanoplot" not in config["skip_stage_dict"]:
        results_list += [*[expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                               datatype=[dat_type, ],
                               stage=["filtered", ],
                               fileprefix=input_file_prefix_dict[dat_type],) for dat_type in long_read_data_type_set],
                        expand(output_dict["qc"] / "nanoplot/{datatype}/{stage}/{fileprefix}.Yield_By_Length.png",
                                   datatype=["nanopore", ],
                                   stage=["trimmed", ],
                                   fileprefix=input_file_prefix_dict["nanopore"],) if "nanopore" in long_read_data_type_set else [],
                        ]

    if config["database_set"]["kraken2"] and kraken_scan_data_type_set and ("kraken" not in config["skip_stage_dict"]):
        results_list += [expand(out_dir_path / "contamination_scan/kraken2/{datatype}/kraken2.{database}.report",
                               datatype=kraken_scan_data_type_set,
                               database=config["database_set"]["kraken2"],
                               )
                        ]



if "filter_draft" in config["stage_list"]:
    results_list += [ ] # TODO: implement
"""
if "draft_qc" in config["stage_list"]:
    stage_dict["contig"] = {}
    assembler_list = ["draft"]
    stage_dict["contig"]["parameters"] = {}
    assembler_option_set_group_dict = {}

    for assembler in assembler_list:
        option_set_group_dict, option_set_group_assignment_dict = None, None
        for option_set in config["coretool_option_sets"][assembler]:
            parameters_label="{0}_{1}".format(assembler, option_set)
            stage_dict["contig"]["parameters"][parameters_label] = {}
            stage_dict["contig"]["parameters"][parameters_label]["included"] = True
            stage_dict["contig"]["parameters"][parameters_label]["assembler"] = assembler
            stage_dict["contig"]["parameters"][parameters_label]["option_set"] = deepcopy(parameters["tool_options"][assembler][option_set])
            if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] is None:
               stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] = config["ploidy"]

            stage_dict["contig"]["parameters"][parameters_label]["haplotype_list"] = ["hap{0}".format(i) for i in range(1, stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] + 1)] if stage_dict["contig"]["parameters"][parameters_label]["option_set"]["assembly_ploidy"] > 1 else ["hap0"]
            stage_dict["contig"]["parameters"][parameters_label]["option_set_group"] = option_set_group_assignment_dict[option_set] if option_set_group_assignment_dict is not None else none

            #for option_supergroup in ["options_affecting_error_correction"]:
            #    stage_dict["contig"]["parameters"][parameters_label][option_supergroup] = option_cluster_reverse_dict[assembler][option_supergroup][option_set]


if "purge_dups" in config["stage_list"]:
    prev_stage = stage_dict["purge_dups"]["prev_stage"]
    purge_dupser_list = config["stage_coretools"]["purge_dups"]["default"]
    stage_dict["purge_dups"]["parameters"] = {}

    for purge_dupser in purge_dupser_list:
        for option_set in config["coretool_option_sets"][purge_dupser]:
            for prev_parameters in stage_dict[prev_stage]["parameters"]:
                parameters_label = "{0}..{1}_{2}".format(prev_parameters, purge_dupser, option_set)
                stage_dict["purge_dups"]["parameters"][parameters_label] = {}
                stage_dict["purge_dups"]["parameters"][parameters_label]["included"] = True
                stage_dict["purge_dups"]["parameters"][parameters_label]["prev_stage"] = prev_stage
                stage_dict["purge_dups"]["parameters"][parameters_label]["prev_parameters"] = prev_parameters
                stage_dict["purge_dups"]["parameters"][parameters_label]["purge_dupser"] = purge_dupser
                stage_dict["purge_dups"]["parameters"][parameters_label]["option_set"] = parameters["tool_options"][purge_dupser][option_set]
                stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"] = stage_dict[stage_dict["purge_dups"]["prev_stage"]]["parameters"][prev_parameters]["haplotype_list"]

    parameters_list = list(stage_dict["purge_dups"]["parameters"].keys())
    results_list += [
                     *[expand(out_dir_path / "purge_dups/{parameters}/{genome_prefix}.purge_dups.{haplotype}.fasta",
                              genome_prefix=[config["genome_prefix"], ],
                              assembly_stage=["contig"],
                              haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                              parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path / "{assembly_stage}/{parameters}/{genome_prefix}.{assembly_stage}.{haplotype}.len",
                             genome_prefix=[config["genome_prefix"], ],
                             assembly_stage=["purge_dups"],
                             haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                             parameters=[parameters_label]) for parameters_label in parameters_list],
                    *[expand(out_dir_path /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/{haplotype}/PB.stat",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["purge_dups"],
                           haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                           parameters=[parameters_label]) for parameters_label in parameters_list],
                    expand(out_dir_path /  "{assembly_stage}/{parameters}/assembly_qc/purge_dups/before.comparison.coverage.png",
                           assembly_stage=["purge_dups"],
                           parameters=parameters_list
                           ),
                    expand(out_dir_path / "{assembly_stage}/{genome_prefix}.{assembly_stage}.stage_stats",
                           genome_prefix=[config["genome_prefix"], ],
                           assembly_stage=["purge_dups"],),
                    [[expand(out_dir_path  / "purge_dups/{parameters}/{purge_stage}/{haplotype}/{genome_prefix}.dups.{artefact}.fasta",
                           purge_stage=["first_stage",] if haplotype == "hap0" else ["first_stage", "second_stage"],
                           genome_prefix=[config["genome_prefix"], ],
                           artefact=["junk", "repeat", "haplotig", "ovlp", "highcov"],
                           haplotype=[haplotype],
                           parameters=[parameters_label]) for haplotype in stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"]] for parameters_label in parameters_list],
                    ]
    if not config["skip_busco"]:
        results_list += [*[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/{genome_prefix}.{assembly_stage}.{haplotype}.busco5.{busco_lineage}.tar.gz",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups", ],
                                haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/haplotype_intersection/{genome_prefix}.{assembly_stage}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                #haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         *[expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/stage_intersection/{genome_prefix}.{haplotype}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                haplotype=stage_dict["purge_dups"]["parameters"][parameters_label]["haplotype_list"],
                                parameters=[parameters_label]) for parameters_label in parameters_list],
                         expand(out_dir_path / "{assembly_stage}/{parameters}/assembly_qc/busco5/all_intersection/{genome_prefix}.{busco_lineage}.busco.merged.tsv",
                                busco_lineage=config["busco_lineage_list"],
                                genome_prefix=[config["genome_prefix"], ],
                                assembly_stage=["purge_dups"],
                                parameters=parameters_list
                                ),
                         ]

if (config["phasing_stage"] in config["stage_list"]) and (not config["skip_phasing"]):

    for datatype in set(data_types) & set(config["read_phasing_data"]):
        if datatype in config["paired_fastq_based_data"]:
            results_list += [*[(expand(out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{assembly_kmer_length}/{datatype}/{pairprefix}_1.fastq.gz",
                                    datatype=[datatype],
                                    stage=[config["phasing_stage"], ],
                                    parameters=[parameters_label],
                                    pairprefix=input_pairprefix_dict[datatype],
                                    genome_prefix=[config["genome_prefix"], ],
                                    haplotype=stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"],
                                    assembly_kmer_length=config["assembly_kmer_length"]
                                    ) if len(stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"]) > 1 else []) for parameters_label in list(stage_dict[config["phasing_stage"]]["parameters"].keys())] ,
                            ]
        else:
            results_list += [*[(expand(out_dir_path / "{stage}/{parameters}/fastq/{haplotype}/{assembly_kmer_length}/{datatype}/{fileprefix}.fastq.gz",
                                    datatype=[datatype],
                                    stage=[config["phasing_stage"], ],
                                    parameters=[parameters_label],
                                    fileprefix=input_file_prefix_dict[datatype],
                                    genome_prefix=[config["genome_prefix"], ],
                                    haplotype=stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"],
                                    assembly_kmer_length=config["assembly_kmer_length"]
                                    ) if len(stage_dict[config["phasing_stage"]]["parameters"][parameters_label]["haplotype_list"]) > 1 else []) for parameters_label in list(stage_dict[config["phasing_stage"]]["parameters"].keys())],
                            ]
"""

#----

#---- Final rule ----
rule all:
    input:
        results_list
        #results_dict[config["mode"]]
#----

#---- Include section ----
include: "workflow/rules/Install/Pip.smk"
include: "workflow/rules/Preprocessing/Files.smk"
include: "workflow/rules/QCFiltering/FastQC.smk"
include: "workflow/rules/QCFiltering/MultiQC.smk"
include: "workflow/rules/QCFiltering/Cutadapt.smk"


#if "nanopore" in data_types:
include: "workflow/rules/QCFiltering/NanoPlot.smk"

include: "workflow/rules/Kmer/Meryl.smk"
include: "workflow/rules/Kmer/Genomescope.smk"


include: "workflow/rules/QCAssembly/BUSCO5.smk"
include: "workflow/rules/QCAssembly/Merqury.smk"
include: "workflow/rules/QCAssembly/QUAST.smk"
include: "workflow/rules/QCAssembly/General.smk"
include: "workflow/rules/Stats/General.smk"
include: "workflow/rules/Contamination/Kraken2.smk"
#include: "workflow/rules/HiC/ReadPhasing.smk"

include: "workflow/rules/Alignment/Index.smk"
include: "workflow/rules/Alignment/Stats.smk"


