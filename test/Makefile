# Test script for CladeCompare
#
# Example foreground and background alignments are taken from:
# Talevich E, Mirza A, Kannan N (2011) Structural and evolutionary divergence of
# eukaryotic protein kinases in Apicomplexa. BMC Evolutionary Biology 11: 321.
# doi:10.1186/1471-2148-11-321.
# https://www.biomedcentral.com/1471-2148/11/321

all: scttlre-v-cdc2.html

clean:
	rm -vf scttlre-v-cdc2.html scttlre-v-cdc2.out scttlre-v-cdc2.pttrn scttlre-v-cdc2-*.pdf

scttlre-v-cdc2.html: ../cladereport.py scttlre-v-cdc2.out
	python $^ -p scttlre-v-cdc2.pttrn > $@

scttlre-v-cdc2.out scttlre-v-cdc2.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	python $^ -p scttlre-v-cdc2.pttrn > scttlre-v-cdc2.out

