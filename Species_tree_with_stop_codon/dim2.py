import re


unknowns = ["Lecanicillium fungicola","Valsa mali","Valsa sordida","Valsa malicola","Valsa mali var. pyri (nom. inval.)","Stachybotrys chlorohalonata IBT 40285"]
sp_list =[]

stps = ["TAG", "TGA", "TAA"]
stp_col ={}
stp_col["TAA"] = "T"
stp_col["TGA"] = "G"
stp_col["TAG"] = "A"

f_in = open("fasta.txt", "r")
seqs = f_in.read()
f_in.close()
seqs = re.sub("Stachybotrys chlorohalonata","Stachybotrys chlorohalonatus",seqs) 

species_stp ={}
seq_list = seqs.split(">")[1:]

prealign = ""
org_list =""

for s in seq_list:
    all_lines = s.splitlines()
    topl = all_lines[0]
    if re.search('"organism_name":"([A-Za-z]+?[^A-Za-z]+?[A-Za-z]+?)[^A-Za-z]', topl):
        fnd_sp = re.search('"organism_name":"([A-Za-z]+?[^A-Za-z]+?[A-Za-z]+?)[^A-Za-z]', topl)
        species = fnd_sp.group(1)
        species1 = re.sub(" ","_",species)
    
        find_oid = re.search('"organism_taxid":"(.*?)"',topl)
        oid = f"{species}\n"
    
    
        seq = "".join(all_lines[1:])
        stp = seq[-3:]
        strt = seq[0:3]
        output = f">{species}:{stp}\n{seq}\n"
        if stp in stps:
            if strt == "ATG":
                if species not in unknowns:
                    if species not in sp_list:
                        prealign = prealign+output
                        org_list = org_list + oid
                        species_stp[species] = stp
                        sp_list.append(species)

f_out = open("Dim2_unaligned.fta", "w")   

f_out.write(prealign)

f_out.close()

f_out = open("Dim2_species.txt", "w")   

f_out.write(org_list[:-1])

f_out.close()     
         
csv_filename = "species_stop.csv"
f_out =open(csv_filename,"w")
f_out.write("species,species_underline,stop,colour\n")
for sp in species_stp:
	stop= species_stp[sp]
	sp1 = re.sub(" ","_",sp)
	f_out.write(f"{sp},{sp1},{species_stp[sp]},{stp_col[stop]}\n")
f_out.close()    
