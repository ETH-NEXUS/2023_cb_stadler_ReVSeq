#!/bin/bash

echo "This script creates the configuration and initial directory structure necessary to run the RevSeq workflow."
echo ""
echo "We need to create and populate six folders in which to store the configuration, the mirror of the metadata files on the sftp server,
        the mirror of the raw data on the sftp server, the raw data in a more user-friendly structure, the anonymized data and the results."
echo "It is strongly suggested this folder is not inside the git repository to avoid to accidentally commit secrets."
echo "The tool will exit if a non-empty folder is provided."
echo ""
src_dir=$(dirname "$0")

read -p "Please provide a the full path where to initialize the configuration (full path): "
    if [[ $REPLY =~ ^[/] ]]
    then
        if [[ -e $REPLY ]]
        then
            echo "Provided path $REPLY exists"
            if [[ ! -d $REPLY ]]
            then
                echo "but it's not a directory. Exiting.."
                exit 200
            fi
            shopt -s nullglob
            shopt -s dotglob 
            # Check for empty files using arrays 
            chk_files=(${REPLY}/*)
            (( ${#chk_files[*]} )) && echo "Directory $REPLY is not empty. Exiting." &&
                    shopt -u nullglob && shopt -u dotglob && exit 201 ||
                    echo "Directory $REPLY is empty."
    	else
            echo "Provided path $REPLY does not exist. Creating the folder"
            mkdir -p $REPLY || exit 500
    	fi
        echo "Initializing configuration"
    mkdir ${REPLY}/config ${REPLY}/raw_data ${REPLY}/metadata ${REPLY}/anondata ${REPLY}/viollier_mirror ${REPLY}/results || exit 301
    cp ${src_dir}/../pipeline_revseq/config_templates/anonymization_table_template.tsv ${REPLY}/config/anonymization_table.tsv || exit 302
    cp ${src_dir}/../pipeline_revseq/config_templates/sample_map_template.tsv ${REPLY}/config/sample_map.tsv || exit 303
    cp ${src_dir}/../pipeline_revseq/config_templates/config_example.yaml ${REPLY}/config/config.yaml || exit 304
    sed -i .bak "s*input_fastqs: \"/path/to/anonymised/fastqs\"*input_fastqs: \"${REPLY}/anondata\"*" ${REPLY}/config/config.yaml
    sed -i .bak "s*sample_map: \"/path/to/sample/map/sample_map.tsv\"*sample_map: \"${REPLY}/config/sample_map.tsv\"*" ${REPLY}/config/config.yaml
    sed -i .bak "s*output_dir: \"path/to/output/dir\"*output_dir: \"${REPLY}/results\"*" ${REPLY}/config/config.yaml
    if [ -f ${REPLY}/config/config.yaml ]
    then
        rm ${REPLY}/config/config.yaml.bak
    fi
    touch ${REPLY}/config/empty_sample.tsv || exit 305
    cp ${src_dir}/../revseq/docker-compose_example.yml ${src_dir}/../revseq/docker-compose.yml || exit 306
    sed -i .bak "s*context: Path/to/revseq/docker*context: $(realpath ${src_dir}/../revseq)*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Raw/Data:/data/raw_data*- ${REPLY}/raw_data:/data/raw_data*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Config:/data/config*- ${REPLY}/config:/data/config*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Results:/data/results*- ${REPLY}/results:/data/results*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Pipeline:/data/pipeline*- $(realpath ${src_dir}/../pipeline_revseq/):/data/pipeline*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Revseq:/data/reveq*- $(realpath ${src_dir}/../revseq):/data/revseq*" ${src_dir}/../revseq/docker-compose.yml
    sed -i .bak "s*- Path/To/Viollier/Mirror:/data/viollier_mirror*- ${REPLY}/viollier_mirror:/data/viollier_mirror*" ${src_dir}/../revseq/docker-compose.yml 
    if [ -f ${src_dir}/../revseq/docker-compose.yml ]
    then
        rm ${src_dir}/../revseq/docker-compose.yml.bak
    fi

    echo "Created the necessary files and directories in ${REPLY}."
    echo "Please follow the rest of the installation instructions to customize the configuration to your specific system"
    exit 0
    else
        echo "The provided string $REPLY does not appear to be an absolute path. Exiting."
        exit 999
    fi
