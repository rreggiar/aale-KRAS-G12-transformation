version 1.0
import '../../../atac.wdl' as atac
import 'compare_md5sum.wdl' as compare_md5sum

workflow test_preseq {
    input {
        File bam
        Boolean paired_end

        File ref_picard_est_lib_size_qc
        File ref_preseq_log

        Float preseq_mem_factor = 0.0
        Float preseq_disk_factor = 5.0
    }

    call atac.preseq { input : 
        paired_end = paired_end,
        bam = bam,
        mem_factor = preseq_mem_factor,
        disk_factor = preseq_disk_factor,
        picard_java_heap = '4G',
    }

    call compare_md5sum.compare_md5sum { input :
        labels = [
            'test_picard_est_lib_size_qc',
            'test_preseq_log',
        ],
        files = select_all([
            preseq.picard_est_lib_size_qc,
            preseq.preseq_log,
        ]),
        ref_files = [
            ref_picard_est_lib_size_qc,
            ref_preseq_log,
        ],
    }
}
