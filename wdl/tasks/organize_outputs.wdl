version 1.0

import "../structs.wdl"

task create_timestamp {
  # Emit a single timestamp shared across the data and index upload passes so both
  # write to the same timestamped bucket directory (keeping each index next to its data).
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
  # Upload one set of outputs (a data pass or an index pass) to the output bucket.
  # upload_outputs.sh writes each file to <bucket>/<workflow>/<version>/<timestamp>/<identifier>/<key>/<basename>,
  # so a data pass and an index pass that share the same keys, identifier, and timestamp land
  # in the same per-key directory.
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

  # Use the base SemVer of the workflow to write outputs
  # e.g. 1.12.0-11-2cdbd238 -> 1.12.0
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
  # Upload one file per sample, each under a key equal to its sample_id, so upload_outputs.sh
  # writes it to <bucket>/<workflow>/<version>/<timestamp>/<identifier>/<sample_id>/<basename>
  # (a per-sample directory). The per-sample-keyed JSON is assembled here from the parallel
  # sample_ids / sample_files arrays because WDL 1.0 cannot build a map with dynamic keys for
  # the caller; upload_outputs.sh itself is used unmodified.
  #
  # sample_files are cloud paths passed as String, NOT File: upload_outputs.sh copies them
  # bucket-to-bucket, so localizing (downloading) them into this task is unnecessary and slow.
  # Declaring them String means Cromwell never localizes them; only their cloud paths are needed.
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

  # Use the base SemVer of the workflow to write outputs
  # e.g. 1.12.0-11-2cdbd238 -> 1.12.0
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
