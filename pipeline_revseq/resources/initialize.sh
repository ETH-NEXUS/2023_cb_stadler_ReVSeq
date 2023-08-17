#!/bin/bash

echo "This script creates the configuration and initial directory structure necessary to run the RevSeq workflow."
echo ""
echo "We need to create and populate six folders in which to store the configuration, the mirror of the metadata files on the sftp server,
        the mirror of the raw data on the sftp server, the raw data in a more user-friendly structure, the anonymized data and the results."
echo "It is strongly suggested this folder is not inside the git repository to avoid to accidentally commit secrets."
echo "The tool will exit if a non-empty folder is provided."
echo ""
src_dir=$(dirname "$0")

read -p "Please provide a the full path where to initialize the configuration (full path): " target_dir
    if [[ $target_dir =~ ^[/] ]]
    then
        if [[ -e $target_dir ]]
        then
            echo "Provided path $target_dir exists"
            if [[ ! -d $target_dir ]]
            then
                echo "but it's not a directory. Exiting.."
                exit 200
            fi
            shopt -s nullglob
            shopt -s dotglob 
            # Check for empty files using arrays 
            chk_files=(${target_dir}/*)
            (( ${#chk_files[*]} )) && echo "Directory $target_dir is not empty. Exiting." &&
                    shopt -u nullglob && shopt -u dotglob && exit 201 ||
                    echo "Directory $target_dir is empty."
    	else
            echo "Provided path $target_dir does not exist. Creating the folder"
            mkdir -p $target_dir || exit 500
    	fi
        echo "Initializing configuration"
    mkdir ${target_dir}/pipeline_configuration ${target_dir}/viollier_mirror ${target_dir}/results ${target_dir}/results/pseudoanon_data || exit 301
    cp ${src_dir}/../../pipeline_revseq/config_template/config_template.yaml ${target_dir}/pipeline_configuration/config.yaml || exit 304
    cp ${src_dir}/../../pipeline_revseq/resources/RespiratoryVirus.20200409.fa ${target_dir}/pipeline_configuration/RespiratoryVirus.20200409.fa ||exit 305

    sed -i "s*input_fastqs: \"/path/to/anonymised/fastqs\"*input_fastqs: \"${target_dir}/anondata\"*" ${target_dir}/pipeline_configuration/config.yaml
    sed -i "s*output_dir: \"path/to/output/dir\"*output_dir: \"${target_dir}/results\"*" ${target_dir}/pipeline_configuration/config.yaml
    sed -i "s*metadata_dir: \"/Path/to/mirror/subdir/with/viollier/metadata\"*output_dir: \"${target_dir}/viollier_mirror/revseq_data\"*" ${target_dir}/pipeline_configuration/config.yaml

    #echo "Extracting the host reference genome. This may take some time."
    #tar -xjf ${src_dir}/../../../config/RespiratoryVirus_hg38.20200409.fa.tar.bz2 --directory ${target_dir}/pipeline_configuration/ || exit 600
    sed -i "s*reference_table: \"/path/to/reference_table\"*reference_table: $(realpath ${src_dir}/../pipeline_revseq/resources/RespiratoryVirus.20200409.bed)\"*" ${target_dir}/pipeline_configuration/config.yaml
    sed -i "s*host_ref: \"/path/to/host/reference.fasta\"*host_ref: $(realpath ${target_dir}/pipeline_configuration/RespiratoryVirus_hg38.20200409.fa)\"*" ${target_dir}/pipeline_configuration/config.yaml
    sed -i "s*reference: \"/path/to/viral/reference/multifasta\"*reference: $(realpath ${target_dir}/pipeline_configuration/RespiratoryVirus.20200409.fa)\"*" ${target_dir}/pipeline_configuration/config.yaml
    sed -i "s*reference_dir: \"/path/to/dir/containing/reference/multifasta\"*reference_dir: $(realpath ${target_dir}/pipeline_configuration/)\"*" ${target_dir}/pipeline_configuration/config.yaml

    cp ${src_dir}/../../revseq/docker-compose_example.yml ${src_dir}/../../revseq/docker-compose.yml || exit 306
    sed -i "s*context: Path/to/revseq/docker*context: $(realpath ${src_dir}/../../revseq)*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Raw/Data:/data/raw_data*- ${target_dir}/raw_data:/data/raw_data*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Config:/data/config*- ${target_dir}/pipeline_configuration:/data/config*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Results:/data/results*- ${target_dir}/results:/data/results*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Pipeline:/data/pipeline*- $(realpath ${src_dir}/../../pipeline_revseq/):/data/pipeline*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Revseq:/data/revseq*- $(realpath ${src_dir}/../revseq):/data/revseq*" ${src_dir}/../../revseq/docker-compose.yml
    sed -i "s*- Path/To/Viollier/Mirror:/data/viollier_mirror*- ${target_dir}/viollier_mirror:/data/viollier_mirror*" ${src_dir}/../../revseq/docker-compose.yml 

    echo ""
    echo "The docker container needs access to the secrets necessary to connect to the Viollier SFTP server."
    read -p "Please provide the location of the ssh config file (full_path): " conf
        if [[ $conf =~ ^[/] ]]
        then
            sed -i "s*file: Path/To/Ssh/config*file: ${conf}*" ${src_dir}/../../revseq/docker-compose.yml
        else
            echo "The provided string ${conf} does not appear to be an absolute path. Exiting."
            exit 999
        fi
    read -p "Please provide the location of the private key file (full_path): " prkey
        if [[ $prkey =~ ^[/] ]]
        then
            sed -i "s*file: Path/To/Private/Key/For/Sftp*file: ${prkey}*" ${src_dir}/../../revseq/docker-compose.yml
        else
            echo "The provided string ${prkey} does not appear to be an absolute path. Exiting."
            exit 999
        fi
    read -p "Please provide the location of the public key file (full_path): " pukey
        if [[ $pukey =~ ^[/] ]]
        then
            sed -i "s*file: Path/To/Public/Key/For/Sftp*file: ${pukey}*" ${src_dir}/../../revseq/docker-compose.yml
        else
            echo "The provided string ${pukey} does not appear to be an absolute path. Exiting."
            exit 999
        fi
    read -p "Please provide the location of the ssh known_hosts file (full_path): " khost
        if [[ $khost =~ ^[/] ]]
        then
        sed -i "s*file: Path/To/Known/Hosts*file: ${khost}*" ${src_dir}/../../revseq/docker-compose.yml
        else
            echo "The provided string ${khost} does not appear to be an absolute path. Exiting."
            exit 999
        fi

    echo ""
    echo "Created the necessary files and directories in ${target_dir}."
    echo "Please follow the rest of the installation instructions to customize the configuration to your specific system"
    exit 0
    else
        echo "The provided string ${target_dir} does not appear to be an absolute path. Exiting."
        exit 999
    fi
