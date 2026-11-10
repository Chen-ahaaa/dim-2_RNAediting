1. Retrieve DIM-2 homologous CDS sequences.
- Download DIM-2 orthologous CDS sequences for Sordariales primarily from OrthoDB.
- Supplement the dataset with additional Sordariales species from JGI, identifying DIM-2 homologs through BLAST searches against their genome sequences.
- Combine the CDS sequences obtained from OrthoDB and JGI into a single FASTA file (fasta.txt).
2. Filter and characterize DIM-2 homologous CDS sequences.
- Use dim2.py to retain one representative sequence per species, assess sequence completeness, and identify stop codon types and positions.
- Generate three output files:
  - Dim2_species.txt: List of non-redundant species.
  - specie_stop.csv: Stop codon information for each species.
  - dim2_unaligned: Filtered, non-redundant DIM-2 CDS sequences in FASTA format.
3. Retrieve the taxonomic tree from NCBI.
- Submit the species list (Dim2_species.txt) to NCBI Taxonomy Common Tree.
- Generate a taxonomy-based tree showing the relationships among the selected Sordariales species.
