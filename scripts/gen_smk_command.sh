sh scripts_for_smk_command/general_cmd_generator.sh | cat metadata/snakemake_prefix_cmd.txt - metadata/snakemake_suffix_cmd.txt | tr '\n' ' ' | awk '{print $0"\n"}'
