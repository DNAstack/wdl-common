version 1.0

import "../structs.wdl"

task trgt {
	input {
		String sample_id
		String? sex

		File bam
		File bam_index

		File reference
		File reference_index
		File tandem_repeat_bed

		String output_prefix

		RuntimeAttributes runtime_attributes
	}

	Boolean sex_defined = defined(sex)
	String karyotype = if select_first([sex, "FEMALE"]) == "MALE" then "XY" else "XX"
	Int threads = 4
	Int disk_size = ceil((size(bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		echo ~{if sex_defined then "" else "Sex is not defined for ~{sample_id}.  Defaulting to karyotype XX for TRGT."}

		trgt --version

		trgt \
			--threads ~{threads} \
			--karyotype ~{karyotype} \
			--genome ~{reference} \
			--repeats ~{tandem_repeat_bed} \
			--reads ~{bam} \
			--output-prefix ~{output_prefix}.trgt

		bcftools --version

		bcftools sort \
			--output-type z \
			--output ~{output_prefix}.trgt.sorted.vcf.gz \
			~{output_prefix}.trgt.vcf.gz

		bcftools index \
			--threads ~{threads - 1} \
			--tbi \
			~{output_prefix}.trgt.sorted.vcf.gz
		
		samtools --version

		samtools sort \
			-@ ~{threads - 1} \
			-o ~{output_prefix}.trgt.spanning.sorted.bam \
			~{output_prefix}.trgt.spanning.bam

		samtools index \
			-@ ~{threads - 1} \
			~{output_prefix}.trgt.spanning.sorted.bam
	>>>

	output {
		File spanning_reads = "~{output_prefix}.trgt.spanning.sorted.bam"
		File spanning_reads_index = "~{output_prefix}.trgt.spanning.sorted.bam.bai"
		File repeat_vcf = "~{output_prefix}.trgt.sorted.vcf.gz"
		File repeat_vcf_index = "~{output_prefix}.trgt.sorted.vcf.gz.tbi"
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/trgt@sha256:88eaa6b6c7d440a48d7f0036e46a2ce4b37cf5be8bd84921eaa69e3c11b98556"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
