# Create Salmon reference
# A. M. Chakrabarti
# revised 13th March 2021

conda activate tdp43-mutants

# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz

grep "^>" <(gunzip -c GRCh38.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

cat gencode.v33.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz

salmon index -t gentrome.fa.gz -d decoys.txt -p 12 -i salmon_index --gencode