version 1.0

import "../structs.wdl"

task hificnv {
	input {
		String sample_id
		String? sex

		File bam
		File bam_index

		File phased_vcf
		File phased_vcf_index

		File reference
		File reference_index

		File exclude_bed
		File exclude_bed_index

		File expected_bed_male
		File expected_bed_female

		String output_prefix

		RuntimeAttributes runtime_attributes
	}

	Boolean sex_defined = defined(sex)
	File expected_bed = if select_first([sex, "FEMALE"]) == "MALE" then expected_bed_male else expected_bed_female

	Int threads = 8
	# Uses ~2 GB memory / thread
	Int mem_gb = threads * 2
	# <1 GB for output
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB"))+ 20)

	command <<<
		set -euo pipefail

		echo ~{if sex_defined then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for HiFiCNV."}

		hificnv --version

		hificnv \
			--threads ~{threads} \
			--bam ~{bam} \
			--ref ~{reference} \
			--maf ~{phased_vcf} \
			--exclude ~{exclude_bed} \
			--expected-cn ~{expected_bed} \
			--output-prefix ~{output_prefix}

		bcftools index --tbi ~{output_prefix}.~{sample_id}.vcf.gz
	>>>

	output {
		File cnv_vcf = "~{output_prefix}.~{sample_id}.vcf.gz"
		File cnv_vcf_index = "~{output_prefix}.~{sample_id}.vcf.gz.tbi"
		File copynum_bedgraph = "~{output_prefix}.~{sample_id}.copynum.bedgraph"
		File depth_bw = "~{output_prefix}.~{sample_id}.depth.bw"
		File maf_bw = "~{output_prefix}.~{sample_id}.maf.bw"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/hificnv@sha256:19fdde99ad2454598ff7d82f27209e96184d9a6bb92dc0485cc7dbe87739b3c2"
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
