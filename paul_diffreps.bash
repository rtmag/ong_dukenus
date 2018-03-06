#
bamToBed -i 1_3502DukeNus_TS543-NT-031117_hg19_i9_rmdup.bam > ../bed/NT_1.bed &
bamToBed -i 2_3502DukeNus_TS543-143-031117_hg19_i10_rmdup.bam > ../bed/sh143_2.bed &
bamToBed -i 3_3502DukeNus_TS543-400-031117_hg19_i11_rmdup.bam > ../bed/sh400_3.bed &
bamToBed -i 4_3502DukeNus_TS543-NT-241117_hg19_i12_rmdup.bam > ../bed/NT_4.bed &
bamToBed -i 5_3502DukeNus_TS543-143-241117_hg19_i13_rmdup.bam > ../bed/sh143_5.bed &
bamToBed -i 6_3502DukeNus_TS543-400-241117_hg19_i14_rmdup.bam > ../bed/sh400_6.bed &

wait

cpan Statistics::TTest
cpan Math::CDF
cpan Parallel::ForkManager

diffReps.pl --treatment ../bed/sh143_2.bed ../bed/sh400_3.bed ../bed/sh143_5.bed ../bed/sh400_6.bed \
--control ../bed/NT_1.bed ../bed/NT_4.bed \
--meth nb --gname hg19 --report ../paul_diffreps/atac.diffreps --frag 0 --nproc 30 

##
##
#
