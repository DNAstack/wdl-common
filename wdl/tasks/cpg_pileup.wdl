version 1.0

import "../structs.wdl"

task cpg_pileup {
	input {
		File bam
		File bam_index

		String output_prefix

		File reference
		File reference_index

		RuntimeAttributes runtime_attributes
	}

	Int threads = 12
	# Uses ~4 GB memory / thread
	Int mem_gb = threads * 4
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		aligned_bam_to_cpg_scores --version

		aligned_bam_to_cpg_scores \
			--threads ~{threads} \
			--bam ~{bam} \
			--ref ~{reference} \
			--output-prefix ~{output_prefix} \
			--min-mapq 1 \
			--min-coverage 10 \
			--model "$PILEUP_MODEL_DIR"/pileup_calling_model.v1.tflite
	>>>

	output {
		Array[File] pileup_beds = glob("~{output_prefix}.*.bed")
		Array[File] pileup_bigwigs = glob("~{output_prefix}.*.bw")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pb-cpg-tools@sha256:b95ff1c53bb16e53b8c24f0feaf625a4663973d80862518578437f44385f509b"
		cpu: threads
		memory: mem_gb + " GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
