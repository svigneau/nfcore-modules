process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high'
    resourceLimits [cpus: 15]

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/9b/9becad054093ad4083a961d12733f2a742e11728fe9aa815d678b882b3ede520/data'
        : 'community.wave.seqera.io/library/cutadapt_trim-galore_pigz:a98edd405b34582d'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    tuple val(meta), path("*report.txt"), emit: log, optional: true
    tuple val(meta), path("*unpaired*.fq.gz"), emit: unpaired, optional: true
    tuple val(meta), path("*.html"), emit: html, optional: true
    tuple val(meta), path("*.zip"), emit: zip, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    // Each --cores N uses: N(read) + N(write) + N(Cutadapt) + 2(extra Cutadapt) + 1(TrimGalore)
    // E.g., --cores 2 uses up to 9 cores, --cores 4 uses up to 15 cores
    def cores = 1
    // Calculate max cores that won't exceed available CPUs
    if (task.cpus) {
        // Divide by ~4 since each core specified multiplies
        cores = (task.cpus as int) / 4
        // Cap at 4 cores due to diminishing returns
        if (cores > 4) {
            cores = 4
        }
        // Ensure at least 1 core
        if (cores < 1) {
            cores = 1
        }
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { it.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads} ${prefix}.fastq.gz
        trim_galore \\
            ${args_list.join(' ')} \\
            --cores ${cores} \\
            --gzip \\
            ${prefix}.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }
    else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trim_galore \\
            ${args} \\
            --cores ${cores} \\
            --paired \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
            cutadapt: \$(cutadapt --version)
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}_trimmed.fq.gz ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.txt"
    }
    else {
        output_command = "echo '' | gzip > ${prefix}_1_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
        output_command += "echo '' | gzip > ${prefix}_2_trimmed.fq.gz ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt"
    }
    """
    ${output_command}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
