#ssu-align
ssu-align SkinMicrobiome_all_seqs.fasta -f SkinMicrobiome_macro-assignment --dna

ssu-mask SkinMicrobiome_macro-assignment   #output SkinMicrobiome_macro-assignment.bacteria.mask.stk

#RAxML
raxmlHPC-PTHREADS-SSE3 -H -f v -m GTRCAT -s SILVA_SkinMicrobiome_TEST_1.bacteria.mask.phy -t LTPs128_SSU_tree_bacteria.newick -T 4 -n LTPs128_SkinMicrobiome

#BLASTn
blastn -query SkinMicrobiome_all_seqs.fasta -out SkinMicrobiome_blast_output_1.txt -task megablast -db 16SMicrobial -outfmt 6  -num_threads 24 -max_target_seqs 1 &


