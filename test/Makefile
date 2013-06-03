# Test script for CladeCompare
#
# Example foreground and background alignments are taken from:
# Talevich E, Mirza A, Kannan N (2011) Structural and evolutionary divergence of
# eukaryotic protein kinases in Apicomplexa. BMC Evolutionary Biology 11: 321.
# doi:10.1186/1471-2148-11-321.
# https://www.biomedcentral.com/1471-2148/11/321

all: scttlre-v-cdc2-gtest.html scttlre-v-cdc2-urn.html scttlre-v-cdc2-jsd.html

clean:
	rm -vf scttlre-v-cdc2-*.html scttlre-v-cdc2-*.out scttlre-v-cdc2-*.pttrn scttlre-v-cdc2-*.pdf

scttlre-v-cdc2-gtest.html: ../cladereport.py scttlre-v-cdc2-gtest.out
	python $^ -p scttlre-v-cdc2-gtest.pttrn > $@

scttlre-v-cdc2-gtest.out scttlre-v-cdc2-gtest.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	python $^ -p scttlre-v-cdc2-gtest.pttrn -s gtest > scttlre-v-cdc2-gtest.out


scttlre-v-cdc2-urn.html: ../cladereport.py scttlre-v-cdc2-urn.out
	python $^ -p scttlre-v-cdc2-urn.pttrn > $@

scttlre-v-cdc2-urn.out scttlre-v-cdc2-urn.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	python $^ -p scttlre-v-cdc2-urn.pttrn -s urn -a .5 > scttlre-v-cdc2-urn.out


scttlre-v-cdc2-jsd.html: ../cladereport.py scttlre-v-cdc2-jsd.out
	python $^ -p scttlre-v-cdc2-jsd.pttrn > $@

scttlre-v-cdc2-jsd.out scttlre-v-cdc2-jsd.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	python $^ -p scttlre-v-cdc2-jsd.pttrn -s jsd -a .5 > scttlre-v-cdc2-jsd.out

