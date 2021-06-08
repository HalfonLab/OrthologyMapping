ORTHOLOGY MAPPINGS FOR TRIBOLIUM, OTHERS?

#================== May 5 2020==========================#

Looks like can get orthology mappings for many species using https://ezmeta.unige.ch/i5k/.
There are a few steps:
(1) X.idmap.txt files map the OrthoDB IDs to species annotation IDs (not sure what happens with species w/o good annotation, if there is such?
Format: ID maps: OrthoDB protein ID to public protein ID mapping ::: SPECIES.idmap.txt

(2) X.97pc.txt files give list of paralogs
Format: longest representative is listed, then comma-separated list of paralogs

(3) A_B.brh files give best reciprocal hits between species A and B
Format: protein_id1 protein_id2 score evalue percent_identity start1 end1 start2 end2


need to convert Tcas OGS3 gene names to OGS2 gene names because orthology mapping only uses the latter.
Did this in two steps:
(1) get all the genes from the 5.2 annotation:
$  perl -ane 'print if $F[2] eq "gene"' GCF_000002335.3_Tcas5.2_genomic.gff > Tcas_5.2_OGS3_genes

(2) then extract any genes with a corresponding 'TC' designation. Other genes will not be included, but they are not included in the orthology mapping file either.

$ perl -ane '$F[8] =~ /ID=(gene\d*);Dbxref=BEETLEBASE:(.*),Gene/; print "$1\t$2\n" unless $seen{$1}++;' Tcas_5.2_OGS3_genes > OGS3toOGS2_conversion.txt


-----------------------
There will need to be several steps to get everything mapped back to Drosophila and then into the SCRMshaw results.
(1) any conversions from current gene set to the gene set used for the orthology mapping
(2) map Species genes to Dmel genes using the appropriate file from i5k/OrthoDB (https://ezmeta.unige.ch/i5k/)
(3) map the Dmel IDs to FlyBase IDs (using DMELA.idmap.txt)
(4) map FlyBase IDs to gene symbols (ftp://ftp.flybase.org/releases/FB2020_02/precomputed_files/synonyms/fb_synonym_fb_2020_02.tsv.gz)

*also should indicate where there a potential paralogs and what their IDs/symbols are