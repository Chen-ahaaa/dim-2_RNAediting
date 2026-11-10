1. extract ortholog list (othodb) for dim2 Sordialles from NCBI -> fasta.txt (https://www.orthodb.org/?ncbi=3876054)
2. parse using dim2.py to extract non-redundant species, check for sequence completeness, find stop -> Dim2_species.txt (species list) and specie_stop.csv, dim2_unaligned.
3. With the species list get the NCBI taxonomy (using common tree) https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi
4. Export to phylip -> phyliptree.phy
5. use NCBI treeviewer to convert to newick format
6. use modify_newick.py to modify file to add in stop codon ->tree_with_stp.nwk
7. Open with TreeView (https://treeviewer.ord
