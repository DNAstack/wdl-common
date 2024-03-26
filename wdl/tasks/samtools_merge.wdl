version 1.0

import "../structs.wdl"

task merge_bams {
	input {
		Array[File] bams

		String output_bam_name

		RuntimeAttributes runtime_attributes
	}

	Int threads = 8
	Int disk_size = ceil(size(bams, "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		samtools --version

		samtools merge \
			-@ ~{threads - 1} \
			-o ~{output_bam_name} \
			~{sep=' ' bams}

		samtools index ~{output_bam_name}
	>>>

	output {
		File merged_bam = "~{output_bam_name}"
		File merged_bam_index = "~{output_bam_name}.bai"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/samtools@sha256:cbe496e16773d4ad6f2eec4bd1b76ff142795d160f9dd418318f7162dcdaa685"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " LOCAL"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
