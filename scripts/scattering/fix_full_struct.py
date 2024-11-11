# xpdb.py -- extensions to Bio.PDB
# (c) 2009 Oliver Beckstein
# Released under the same license as Biopython.
# See http://biopython.org/wiki/Reading_large_PDB_files
#%%
import sys
import Bio.PDB
import Bio.PDB.StructureBuilder
def main(infname, outfname):
  # Serial numbers are taken as remainders of 1e4
  folder = "targets/"
  with open(folder+infname,'r') as infile, open(folder+outfname,'w') as outfile:
    line = infile.readline()
    while line!="":
        if len(line) < 18: 
           line = infile.readline()
           continue
        while line[18] == " ": # Incorrect space, need to remove.
            line = line[0:18] + line[19:]
        outfile.write(line)
        line = infile.readline()


main("4et8_full_struct.pdb","4et8_full_struct_fixed.pdb")
# %%
