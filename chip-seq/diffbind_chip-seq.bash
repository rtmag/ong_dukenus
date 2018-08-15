macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/shNT-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n shNT_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/shNT_mnase_2_rmdup.bam &

#
macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh143_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/sh143_input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh143_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/sh143_mnase_2_rmdup.bam &

#
macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh400_1 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam -c /root/ong_dukenus/chip-seq/bam/sh400-input_1_rmdup.bam &

macs2 callpeak -f BAMPE -g hs -q 0.00001 --broad --keep-dup auto -n sh400_2 --outdir /root/ong_dukenus/chip-DIFF/macs/ \
-t /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam -c /root/ong_dukenus/mnase_batch2/bam/sh400_mnase_2_rmdup.bam &
#
# DIFFREPS
bamToBed -i /root/ong_dukenus/chip-seq/bam/shNT-IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/shNT_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/shNT_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/shNT_IP_2_rmdup.bed &
bamToBed -i /root/ong_dukenus/chip-seq/bam/sh143_IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh143_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/sh143_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh143_IP_2_rmdup.bed &
bamToBed -i /root/ong_dukenus/chip-seq/bam/sh400-IP_1_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh400_IP_1_rmdup.bed &
bamToBed -i /root/ong_dukenus/mnase_batch2/bam/sh400_IP_2_rmdup.bam > /root/ong_dukenus/chip-DIFF/bed/sh400_IP_2_rmdup.bed &


diffReps.pl --treatment ./bed/sh143_IP_1_rmdup.bed ./bed/sh143_IP_2_rmdup.bed \
./bed/sh400_IP_1_rmdup.bed ./bed/sh400_IP_2_rmdup.bed \
--control ./bed/shNT_IP_1_rmdup.bed ./bed/shNT_IP_2_rmdup.bed \
--meth nb --gname hg19 --report /root/ong_dukenus/chip-DIFF/diffreps/diffChip_both --frag 0 --nproc 30 &

diffReps.pl --treatment ./bed/sh143_IP_2_rmdup.bed \
./bed/sh400_IP_2_rmdup.bed \
--control ./bed/shNT_IP_2_rmdup.bed \
--meth gt --gname hg19 --report /root/ong_dukenus/chip-DIFF/diffreps/diffChip_batch2 --frag 0 --nproc 30 &
