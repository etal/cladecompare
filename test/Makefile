# Test script for CladeCompare
#
# Example foreground and background alignments are taken from:
# Talevich E, Mirza A, Kannan N (2011) Structural and evolutionary divergence of
# eukaryotic protein kinases in Apicomplexa. BMC Evolutionary Biology 11: 321.
# doi:10.1186/1471-2148-11-321.
# https://www.biomedcentral.com/1471-2148/11/321

PY=python

all: scttlre-v-cdc2-gtest.html scttlre-v-cdc2-urn.html scttlre-v-cdc2-jsd.html \
	scttlre-v-cdc2-hypg.html cdc2-gtest.html cdc2-jsd.html cdc2-urn.html \
	cdc2-hypg.html 

clean:
	rm -vf scttlre-v-cdc2-*.html scttlre-v-cdc2-*.out scttlre-v-cdc2-*.pttrn \
		scttlre-v-cdc2-*.pdf cdc2-*.html cdc2-*.out cdc2-*.pttrn


# Paired comparisons

scttlre-v-cdc2-gtest.html: ../cladereport.py scttlre-v-cdc2-gtest.out
	$(PY) $^ -p scttlre-v-cdc2-gtest.pttrn -o $@

scttlre-v-cdc2-gtest.out scttlre-v-cdc2-gtest.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	$(PY) $^ -p scttlre-v-cdc2-gtest.pttrn -s gtest -w -o scttlre-v-cdc2-gtest.out


scttlre-v-cdc2-urn.html: ../cladereport.py scttlre-v-cdc2-urn.out
	$(PY) $^ -p scttlre-v-cdc2-urn.pttrn -o $@

scttlre-v-cdc2-urn.out scttlre-v-cdc2-urn.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	$(PY) $^ -p scttlre-v-cdc2-urn.pttrn -s urn -o scttlre-v-cdc2-urn.out


scttlre-v-cdc2-jsd.html: ../cladereport.py scttlre-v-cdc2-jsd.out
	$(PY) $^ -p scttlre-v-cdc2-jsd.pttrn -o $@

scttlre-v-cdc2-jsd.out scttlre-v-cdc2-jsd.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	$(PY) $^ -p scttlre-v-cdc2-jsd.pttrn -s jsd -a .2 -o scttlre-v-cdc2-jsd.out

scttlre-v-cdc2-hypg.html: ../cladereport.py scttlre-v-cdc2-hypg.out
	$(PY) $^ -p scttlre-v-cdc2-hypg.pttrn -o $@

scttlre-v-cdc2-hypg.out scttlre-v-cdc2-hypg.pttrn: ../cladecompare.py scttlre-domain.fasta_aln.cma cdc2-domain.fasta_aln.cma
	$(PY) $^ -p scttlre-v-cdc2-hypg.pttrn -s hypg -o scttlre-v-cdc2-hypg.out


# Single mode

cdc2-gtest.html: ../cladereport.py cdc2-gtest.out
	$(PY) $^ -p cdc2-gtest.pttrn -o $@

cdc2-gtest.out cdc2-gtest.pttrn: ../cladecompare.py cdc2-domain.fasta_aln.cma
	$(PY) $^ -p cdc2-gtest.pttrn -s gtest -w -a .001 -o cdc2-gtest.out


cdc2-urn.html: ../cladereport.py cdc2-urn.out
	$(PY) $^ -p cdc2-urn.pttrn -o $@

cdc2-urn.out cdc2-urn.pttrn: ../cladecompare.py cdc2-domain.fasta_aln.cma
	$(PY) $^ -p cdc2-urn.pttrn -s urn -o cdc2-urn.out


cdc2-jsd.html: ../cladereport.py cdc2-jsd.out
	$(PY) $^ -p cdc2-jsd.pttrn -o $@

cdc2-jsd.out cdc2-jsd.pttrn: ../cladecompare.py cdc2-domain.fasta_aln.cma
	$(PY) $^ -p cdc2-jsd.pttrn -s jsd -a .2 -o cdc2-jsd.out

cdc2-hypg.html: ../cladereport.py cdc2-hypg.out
	$(PY) $^ -p cdc2-hypg.pttrn -o $@

cdc2-hypg.out cdc2-hypg.pttrn: ../cladecompare.py cdc2-domain.fasta_aln.cma
	$(PY) $^ -p cdc2-hypg.pttrn -s hypg -a .2 -o cdc2-hypg.out

