#script_sampling_HCC
cut -d ' ' -f1 BVBRC_genome_sequence_subset.fasta | sed 's/|/_/g' > new_subset.fa
conda activate barrnap
mkdir individual_genomes
mv accn_* individual_genomes
cd individual_genomes
for i in *.fasta; do echo $i; barrnap $i --outseq barrnap_$i; done
cat barrnap_* > all_samples.fasta
grep ">" all_samples.fasta | sed 's/>//g' > headers_barrnap.txt
grep "16S" headers_barrnap.txt > headers_16S.txt
conda activate seqkit
for i in barrnap_*; do echo $i; seqkit grep -n -f headers_16S.txt $i > 16S_$i; done

mkdir 16S_files
mv 16S_*.fasta 16S_files
cd 16S_files
find . -name '*.fasta' -size 0 -print0 | xargs -0 rm


conda activate badread
mkdir simulated_reads_human_default
for i in *.fasta; do echo $i; badread simulate --reference $i --quantity 5M --length 1500,200 --qscore_model nanopore2020 \ 
--glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --start_adapter_seq "" --end_adapter_seq "" > simulated_reads_human_default/sim_$i.fastq; done

cd simulated_reads_human_default

conda activate nanofilt
mkdir filtered
for i in *.fastq; do echo $i; NanoFilt -q 10 -l 1000 --maxlength 1700 $i > filtered/filt_$i; done

cd filtered

conda activate fastq-tools
mkdir sampled
for i in *.fastq; do echo $i; fastq-sample -n 2000  -r -s 45 -o sampled/samp_$i $i; done

cd sampled/
find . -name '*fastq' -exec bash -c ' mv $0 ${0/\.fasta.fastq.fastq/.fastq}' {} \;

