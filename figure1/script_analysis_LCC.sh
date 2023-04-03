

#download ZymoBiomics community standard fasta files from the Zymo web and decompress
##simulate ONT reads using badread
conda activate badread
mkdir ident_default
for i in *.fasta; do echo $i; badread simulate --reference $i --quantity 40M --length 1500,200 --qscore_model nanopore2020 --glitches 0,0,0 --junk_reads 0 --random_reads 0 --chimeras 0 --start_adapter_seq "" \
--end_adapter_seq "" > ident_default/sim_$i.fastq; done

cd ident_default
conda activate nanofilt
for i in *.fastq; do echo $i; Nanofilt -q 10 -l 1000 --maxlength 1700 $i > filt_$i; done

mkdir simulate
for i in *Escherichia_coli_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 4040 -r -s $f -o simulate/$f$i $i; done; done
for i in *Bacillus_subtilis_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 6960 -r -s $f -o simulate/$f$i $i; done; done
for i in *Enterococcus_faecalis_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 3960 -r -s $f -o simulate/$f$i $i; done; done
for i in *Lactobacillus_fermentum_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 7360 -r -s $f -o simulate/$f$i $i; done; done
for i in *Pseudomonas_aeruginosa_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 1680 -r -s $f -o simulate/$f$i $i; done; done
for i in *Salmonella_enterica_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 4160 -r -s $f -o simulate/$f$i $i; done; done
for i in *Staphylococcus_aureus_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 6200 -r -s $f -o simulate/$f$i $i; done; done
for i in *Listeria_monocytogenes_*; do echo $i; for f in {1..10}; do echo "$f"; fastq-sample -n 5640 -r -s $f -o simulate/$f$i $i; done; done

cd simulate
mkdir mock{1..10}
mv 10* mock10/
mv 1* mock1/
mv 2* mock2/
mv 3* mock3/
mv 4* mock4/
mv 5* mock5/
mv 6* mock6/
mv 7* mock7/
mv 8* mock8/
mv 9* mock9/

for i in mock*
    do
    echo $i
    cat $i/*.fastq > $i/$i.fastq
done

mkdir all_mocks
cp mock*/mock*.fastq all_mocks


##########################
####wf-metagenomics#######
##########################

cd all_mocks
#MANUAL STEP: create a file tree for the simulated reads that resemble a ONT run
#in the fastq_pass directory ten directories barcode[0-10] were created and simulated reads were move to each corresponding number (mock1.fastq was moved to barcode01 and so on)
conda activate nextflow
nextflow run epi2me-labs/wf-metagenomics --fastq fastq_pass/ --classifier minimap2 --min_len 1000 --max_len 1700 --threads 8 -profile singularit

##########################
########EMU###############
##########################

#in the case of EMU single mock files were used
for i in mock*; do echo $i; emu abundance --db /mnt/cive/csalazar/databases/emu_silva/ --type map-ont --output-dir emu_output --keep-counts --output-unclassified $i  --threads 32; done


##########################
#####porefile#############
##########################
#for porefile single mock files were used

nextflow run microgenlab/porefile -r polishing -with-singularity /mnt/cive/csalazar/software/porefile_poneleunnombre.sif --fq 'mock*.fastq' --outdir results \
--nanofilt_quality 10 --nanofilt_length 1000 --nanofilt_maxlength 1700 --megan_lcaAlgorithm naive -profile nagual --isDemultiplexed --lowAbundanceThreshold 0.025 --megan_topPercentPolish 5

mkdir raw
cp -r output emu_output results

#run script_figure1.R to analze and generate visualizations
