version 1.0

task star {
    input {
        String SID 
        File star_index
        File fastq1
        File fastq2
        # STAR options
        String outFilterType = "BySJout"
        String outSAMtype = "BAM SortedByCoordinate"
        String quantMode = "TranscriptomeSAM"
        String outSAMattributes = "NH HI AS NM MD nM"
        
        Int memory
        Int disk_space
        Int ncpu

        String docker
    }

    command <<<
        
        echo "FASTQs:"
        echo "$FASTQ_1_ABS"
        echo "$FASTQ_2_ABS"

        # extract index
        echo "--- $(date "+[%b %d %H:%M:%S]") Extracting STAR index ---"
        mkdir star_index
        tar -xvf ~{star_index} -C star_index --strip-components=1
        echo "--- $(date "+[%b %d %H:%M:%S]")" Done extracting index

        mkdir star_out
        echo "--- $(date "+[%b %d %H:%M:%S]") Running STAR ---"
        STAR --genomeDir star_index \
            --readFilesIn "$FASTQ_1_ABS" "$FASTQ_2_ABS" \
            --outFileNamePrefix star_out/~{SID}. \
            --readFilesCommand zcat \
            --outSAMattributes ~{outSAMattributes} \
            --outFilterType ~{outFilterType} \
            --runThreadN ~{ncpu} \
            --outSAMtype ~{outSAMtype} \
            --quantMode ~{quantMode} \
            --outSAMunmapped Within \
            --outSAMstrandField intronMotif \
            --outSAMattrRGline ID:rg1 PL:Illumina LB:~{SID} SM:~{SID}

        cd star_out
        ls

        echo "--- $(date "+[%b %d %H:%M:%S]") Running samtools index ---"
        samtools index ~{SID}.Aligned.sortedByCoord.out.bam

        echo "--- $(date "+[%b %d %H:%M:%S]") Finished running samtools index ---"
        ls
        echo "--- $(date "+[%b %d %H:%M:%S]") Completed task ---"
    >>>

    output {
        File bam_file = "star_out/${SID}.Aligned.sortedByCoord.out.bam"
        File bam_index = "star_out/${SID}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "star_out/${SID}.Aligned.toTranscriptome.out.bam"
        File junctions = "star_out/${SID}.SJ.out.tab"
        Array[File] logs = [
            "star_out/${SID}.Log.final.out",
            "star_out/${SID}.Log.out",
            "star_out/${SID}.Log.progress.out"
        ]
    }

    runtime {
        singularity: "docker://${docker}"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${ncpu}"

    }

    parameter_meta {
        SID: {
            type: "id"
        }
        fastq1: {
            label: "Forward End Read FASTQ File"
        }
        fastq2: {
            label: "Reverse End Read FASTQ File"
        }
        star_index: {
            label: "Genome STAR Index File"
        }
    }

    meta {
        author: "Archana Raja"
    }
}
