version 1.0

import "../structs.wdl"

task paraphase {
	input {
		File bam
		File bam_index

		File reference
		File reference_index

		String sample_id
		String out_directory

		RuntimeAttributes runtime_attributes
	}

	Int threads = 8
	Int mem_gb = 16
	Int disk_size = ceil(size(bam, "GB") + 20)

	command <<<
		set -euo pipefail

		paraphase --version

		paraphase \
			--threads ~{threads} \
			--bam ~{bam} \
			--reference ~{reference} \
			--out ~{out_directory}

		if ls ~{out_directory}/~{sample_id}_vcfs/*.vcf &> /dev/null; then
			cd ~{out_directory} \
				&& tar zcvf ~{out_directory}.tar.gz ~{sample_id}_vcfs/*.vcf \
				&& mv ~{out_directory}.tar.gz ../
		fi
	>>>

	output {
		File output_json = "~{out_directory}/~{sample_id}.json"
		File realigned_bam = "~{out_directory}/~{sample_id}_realigned_tagged.bam"
		File realigned_bam_index = "~{out_directory}/~{sample_id}_realigned_tagged.bam.bai"
		File? paraphase_vcfs = "~{out_directory}.tar.gz"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/paraphase@sha256:b9852d1a43485b13c563aaddcb32bacc7f0c9088c2ca007051b9888e9fe5617d"
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