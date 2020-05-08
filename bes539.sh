#!/bin/bash
#SBATCH -t 120:00:00
#SBATCH -J Single_copy_orthologs
#SBATCH --mem=120G
#SBATCH --exclusive

module load SRA-Toolkit/2.9.0-centos_linux64
module load FastQC/0.11.8-Java-1.8 
module load Trimmomatic/0.38-Java-1.8 
module load SPAdes/3.11.1-foss-2016b 
module load TransDecoder/5.5.0 
module load hmmer/3.1b
module load MCL/14.137-foss-2018b 
module load OrthoFinder/2.3.3-foss-2018b-Python-2.7.15
module load DIAMOND/0.9.25-foss-2018b

cd /data3/lanelab/terpis/BES539
### Get Data from NCBI SRA
fastq-dump --split-files SRR1294384 #Dinobryon
fastq-dump --split-files SRR1300417 #cafeteria
fastq-dump --split-files SRR1300413 #aplanochytrium
fastq-dump --split-files SRR1294406 #Cylindrotheca

### Rename
mv SRR1294384_1* Dinobryon_1.fq 
mv SRR1294384_2* Dinobryon_2.fq
mv SRR1300417_1* Cafe_1.fq
mv SRR1300417_2* Cafe_2.fq
mv SRR1300413_1* Aplano_1.fq
mv SRR1300413_2* Aplano_2.fq
mv SRR1294406_1* Cylindro_1.fq
mv SRR1294406_2* Cylindro_2.fq

ls *.fq > stramenopiles.txt

### QualityReport
for f in $(cat stramenopiles.txt); do fastqc $f; done

### Trimming
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 40 \
-validatePairs \
Dinobryon_1.fq \
Dinobryon_2.fq \
Dinobryon_forward_paired.fq Dinobryon_forward_unpaired.fq \
Dinobryon_reverse_paired.fq Dinobryon_reverse_unpaired.fq \
ILLUMINACLIP:/data3/lanelab/terpis/BES539/adapters/TruSeq3-PE-2.fa:2:30:10 \
HEADCROP:15 \
LEADING:30 \
TRAILING:30 \
MINLEN:30

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 40 \
-validatePairs \
Cylindro_1.fq \
Cylindro_2.fq \
Cylindro_forward_paired.fq Cylindro_forward_unpaired.fq \
Cylindro_reverse_paired.fq Cylindro_reverse_unpaired.fq \
ILLUMINACLIP:/data3/lanelab/terpis/BES539/adapters/TruSeq3-PE-2.fa:2:30:10 \
HEADCROP:15 \
LEADING:30 \
TRAILING:30 \
MINLEN:80

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 40 \
-validatePairs \
Aplano_1.fq \
Aplano_2.fq \
Aplano_forward_paired.fq Aplano_forward_unpaired.fq \
Aplano_reverse_paired.fq Aplano_reverse_unpaired.fq \
ILLUMINACLIP:/data3/lanelab/terpis/BES539/adapters/TruSeq3-PE-2.fa:2:30:10 \
HEADCROP:15 \
LEADING:30 \
TRAILING:30 \
MINLEN:30

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads 40 \
-validatePairs \
Cafe_1.fq \
Cafe_2.fq \
Cafe_forward_paired.fq Cafe_forward_unpaired.fq \
Cafe_reverse_paired.fq Cafe_reverse_unpaired.fq \
ILLUMINACLIP:/data3/lanelab/terpis/BES539/adapters/TruSeq3-PE-2.fa:2:30:10 \
HEADCROP:15 \
LEADING:30 \
TRAILING:30 \
MINLEN:30
### Assembly
rnaspades.py --pe1-1 Cafe_forward_paired.fq --pe1-2 Cafe_reverse_paired.fq -t 60 -o Cafe_rnaSPAdes
rnaspades.py --pe1-1 Aplano_forward_paired.fq --pe1-2 Aplano_reverse_paired.fq -t 60 -o Cafe_rnaSPAdes
rnaspades.py --pe1-1 Cylindro_forward_paired.fq --pe1-2 Cylindro_reverse_paired.fq -t 60 -o Cafe_rnaSPAdes
rnaspades.py --pe1-1 Dinobryon_forward_paired.fq --pe1-2 Dinobryon_reverse_paired.fq -t 60 -o Cafe_rnaSPAdes

mkdir Assembled_transcripts

mv Aplano_rnaSPAdes/transcripts.fasta Assembled_transcripts/Aplano.fasta
mv Cafe_rnaSPAdes/transcripts.fasta Assembled_transcripts/Cafe.fasta
mv Cylindro_rnaSPAdes/transcripts.fasta Assembled_transcripts/Cylindro.fasta
mv Dinobryon_rnaSPAdes/transcripts.fasta Assembled_transcripts/Dinobryon.fasta

cd Assembled_transcripts
ls | grep "fasta$" | cut -f1 -d"." > data.txt #this creates a txt file containing only species names

###Translate
for f in $(cat data.txt); do TransDecoder.LongOrfs -t "$f".fasta; done
for f in $(cat data.txt); do TransDecoder.Predict -t "$f".fasta; done

###Renaming
for f in *.pep; do sed -i "s/ /_/g; s/~~//g; s/,//g; s/://g; s/-//g; s/#//g; s/(//g; s/)//g; s/\[/___/g; s/\]//g; s/\./_/g; s/\//_/g; s/;//g; s/\'//g; s/@/_/g; s/=//g" $f; done #remove any illegal characters
for f in *.pep; do awk '{if($0 ~ /^>/){x= substr($0, 1, 200); print x} else {print $0}}' $f > $f.short; done  #limits length of header to 200 characters
sed 's/NODE.*_g/Cafe_g/' Cafe.fasta.transdecoder.pep.short > Cafe.faa
sed 's/NODE.*_g/Aplano_g/' Aplano.fasta.transdecoder.pep.short > Aplano.faa
sed 's/NODE.*_g/Dinobryon_g/' Dinobryon.fasta.transdecoder.pep.short > Dinobryon.faa
sed 's/NODE*_g/Cylindro_g/' Cylindro.fasta.transdecoder.pep.short > Cylindro.faa

###Single Copy orthologs
python /opt/software/OrthoFinder/2.3.3-foss-2018b-Python-2.7.15/orthofinder.py -f /data3/lanelab/terpis/BES539/Assembled_transcripts/OrthoFinder/

