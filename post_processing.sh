
java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/1_3502DukeNus_TS543-NT-031117_hs_i9_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/1_3502DukeNus_TS543-NT-031117_hs_i9_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/1_3502DukeNus_TS543-NT-031117_hs_i9_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/3_3502DukeNus_TS543-400-031117_hs_i11_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/3_3502DukeNus_TS543-400-031117_hs_i11_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/3_3502DukeNus_TS543-400-031117_hs_i11_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/5_3502DukeNus_TS543-143-241117_hs_i13_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/5_3502DukeNus_TS543-143-241117_hs_i13_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/5_3502DukeNus_TS543-143-241117_hs_i13_Aligned.sortedByCoord.out.mfile

java -jar /root/myPrograms/picard/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
I=/root/ong_dukenus/bam/6_3502DukeNus_TS543-400-241117_hs_i14_Aligned.sortedByCoord.out.bam \
O=/root/ong_dukenus/bam/6_3502DukeNus_TS543-400-241117_hs_i14_Aligned.sortedByCoord.rmdup.out.bam \
M=/root/ong_dukenus/bam/6_3502DukeNus_TS543-400-241117_hs_i14_Aligned.sortedByCoord.out.mfile

samtools index /root/ong_dukenus/bam/1_3502DukeNus_TS543-NT-031117_hs_i9_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/3_3502DukeNus_TS543-400-031117_hs_i11_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/4_3502DukeNus_TS543-NT-241117_hs_i12_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/5_3502DukeNus_TS543-143-241117_hs_i13_Aligned.sortedByCoord.rmdup.out.bam &
samtools index /root/ong_dukenus/bam/6_3502DukeNus_TS543-400-241117_hs_i14_Aligned.sortedByCoord.rmdup.out.bam &


