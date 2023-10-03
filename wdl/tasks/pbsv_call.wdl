version 1.0

# Call SVs using pbsv

import "../structs.wdl"

task pbsv_call {
	input {
		String sample_id
		Array[File] svsigs
		Int? sample_count

		File reference
		File reference_index
		String reference_name

		Array[String] regions

		Int mem_gb = if select_first([sample_count, 1]) > 3 then 96 else 64

		RuntimeAttributes runtime_attributes
	}

	Int threads = 8
	Int disk_size = ceil((size(svsigs, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		# pbsv has the ability to call SVs by region by using indexed signatures, but
		#   if an svsig.gz file doesn't contain any signatures in the region, then
		#   pbsv crashes. To avoid this, filter the svsig.gz files to only contain
		#   signatures in the regions.
		# This is brittle and likely to break if pbsv discover changes output format.
		# Build a pattern to match; we want headers (e.g., '^#') and signature
		#   records where third column matches the chromosome (e.g., '^.\t.\tchr1\t').
		pattern=$(echo ~{sep=" " regions} \
			| sed 's/^/^.\\t.\\t/; s/ /\\t\|^.\\t.\\t/; s/$/\\t/' \
			| echo "^#|""$(</dev/stdin)")

		for svsig in ~{sep=" " svsigs}; do
			svsig_basename=$(basename "$svsig" .svsig.gz)
			gunzip -c "$svsig" \
				| grep -P "$pattern" \
				| bgzip -c > "${svsig_basename}.regions.svsig.gz" \
				&& echo "${svsig_basename}.regions.svsig.gz" >> svsigs.fofn
		done

		pbsv --version

		pbsv call \
			--hifi \
			--min-sv-length 20 \
			--log-level INFO \
			--num-threads ~{threads} \
			~{reference} \
			svsigs.fofn \
			"~{sample_id}.~{reference_name}.pbsv.vcf"
	>>>

	output {
		Array[File] pbsv_vcf = glob("~{sample_id}.~{reference_name}.pbsv.vcf")
	}

	runtime {
		docker: "~{runtime_attributes.container_registry}/pbsv@sha256:d78ee6deb92949bdfde98d3e48dab1d871c177d48d8c87c73d12c45bdda43446"
		cpu: threads
		memory: "~{mem_gb} GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: runtime_attributes.preemptible_tries
		maxRetries: runtime_attributes.max_retries
		awsBatchRetryAttempts: runtime_attributes.max_retries
		queueArn: runtime_attributes.queue_arn
		zones: runtime_attributes.zones
	}
}
