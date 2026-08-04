version 1.0

import "../structs.wdl"

task create_timestamp {
  input {
    RuntimeAttributes runtime_attributes
  }

  String tools_docker_image = (if (runtime_attributes.container_registry == "quay.io/pacbio") then "dnastack/" else runtime_attributes.container_registry + "/") + "hifi_solves_tools:2.1.2"

  command <<<
    set -euo pipefail

    date +%s > timestamp.txt
  >>>

  output {
    String timestamp = read_string("timestamp.txt")
  }

  runtime {
    docker: tools_docker_image
    cpu: 2
    memory: "4 GB"
    disks: "local-disk 15 HDD"
    disk: "15 GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpu_platform
  }
}

task organize_outputs {
  input {
    String identifier

    Map[String, Array[String]] workflow_outputs_pre_upload

    String timestamp

    String workflow_name
    String workflow_version

    String output_bucket
    String backend

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  String workflow_version_semver = sub(workflow_version, "-.*$", "")
  String tools_docker_image = (if (runtime_attributes.container_registry == "quay.io/pacbio") then "dnastack/" else runtime_attributes.container_registry + "/") + "hifi_solves_tools:2.1.2"

  command <<<
    set -euo pipefail

    < ~{write_json(workflow_outputs_pre_upload)} \
      jq 'if type == "array" then map({(.left): .right}) | add else . end' \
      > "~{identifier}.outputs.json"

    upload_outputs.sh \
      -b "~{backend}" \
      -i "~{identifier}" \
      -t "~{timestamp}" \
      -w "~{workflow_name}" \
      -v "~{workflow_version_semver}" \
      -o "~{output_bucket}" \
      -j "~{identifier}.outputs.json" \
      -p "~{identifier}.outputs"
  >>>

  output {
    File output_manifest_json = "~{identifier}.outputs.manifest.json"
  }

  runtime {
    docker: tools_docker_image
    cpu: 2
    memory: "4 GB"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpu_platform
  }
}

task organize_per_sample_outputs {
  input {
    String identifier

    Array[String] sample_ids
    Array[String] sample_files

    String timestamp

    String workflow_name
    String workflow_version

    String output_bucket
    String backend

    RuntimeAttributes runtime_attributes
  }

  Int disk_size = 20

  String workflow_version_semver = sub(workflow_version, "-.*$", "")
  String tools_docker_image = (if (runtime_attributes.container_registry == "quay.io/pacbio") then "dnastack/" else runtime_attributes.container_registry + "/") + "hifi_solves_tools:2.1.2"

  command <<<
    set -euo pipefail

    # Build {"<sample_id>": ["<file>"], ...} from the parallel sample_ids / sample_files arrays.
    paste -d '\t' ~{write_lines(sample_ids)} ~{write_lines(sample_files)} > sample_outputs.tsv
    echo '{}' > "~{identifier}.outputs.json"
    while IFS="$(printf '\t')" read -r sample_id sample_file; do
      jq --arg key "${sample_id}" --arg path "${sample_file}" \
        '. + {($key): [$path]}' "~{identifier}.outputs.json" > "~{identifier}.outputs.json.tmp"
      mv "~{identifier}.outputs.json.tmp" "~{identifier}.outputs.json"
    done < sample_outputs.tsv

    upload_outputs.sh \
      -b "~{backend}" \
      -i "~{identifier}" \
      -t "~{timestamp}" \
      -w "~{workflow_name}" \
      -v "~{workflow_version_semver}" \
      -o "~{output_bucket}" \
      -j "~{identifier}.outputs.json" \
      -p "~{identifier}.outputs"
  >>>

  output {
    File output_manifest_json = "~{identifier}.outputs.manifest.json"
  }

  runtime {
    docker: tools_docker_image
    cpu: 2
    memory: "4 GB"
    disks: "local-disk ~{disk_size} HDD"
    disk: "~{disk_size} GB"
    preemptible: runtime_attributes.preemptible_tries
    maxRetries: runtime_attributes.max_retries
    zones: runtime_attributes.zones
    cpuPlatform: runtime_attributes.cpu_platform
  }
}
