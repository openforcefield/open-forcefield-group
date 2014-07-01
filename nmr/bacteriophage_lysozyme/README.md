======================
Bacteriophage Lysozyme


http://www.bmrb.wisc.edu/data_library/summary/index.php?bmrbId=16664


Note: 1AM7 has the unnatural amino acid NZ2-TRP (TRN), which I manually edited to be TRP


NOTES: Because of the low pH, we use the following protocol to generate the starting structure:

1.  Pass 1AM7 through PDBFixer to fix TRN-TRP issue
2.  Use H++ to guess the pKa of all residues and assign to maximum probability state at pH 5.45
3.  Use output PDB as input to pipeline.

This leads to 2 / 3 HIS residues protonated, but the buried HIS48 is in the neutral state.



http://rcsb.org/pdb/explore/explore.do?structureId=1AM7



      'ionic strength'   0     .   M   16664 1 
       pH                5.45 0.05 pH  16664 1 
       pressure          1     .   atm 16664 1 
       temperature     293    0.5  K   16664 1 


>1AM7:A|PDBID|CHAIN|SEQUENCE
MVEINNQRKAFLDMLAWSEGTDNGRQKTRNHGYDVIVGGELFTDYSDHPRKLVTLNPKLKSTGAGRYQLLSRWWDAYRKQ
LGLKDFSPKSQDAVALQQIKERGALPMIDRGDIRQAIDRCSNIWASLPGAGYGQFEHKADSLIAKFKEAGGTVREIDV


16664 sequence:

MVEINNQRKAFLDMLAWSEG
TDNGRQKTRNHGYDVIVGGE
LFTDYSDHPRKLVTLNPKLK
STGAGRYQLLSRWWDAYRKQ
LGLKDFSPKSQDAVALQQIK
ERGALPMIDRGDIRQAIDRC
SNIWASLPGAGYGQFEHKAD
SLIAKFKEAGGTVREIDV


19127 sequence:
MVEINNQRKAFLDMLAWSEG
TDNGRQKTRNHGYDVIVGGE
LFTDYSDHPRKLVTLNPKLK
STGAGRYQLLSRWWDAYRKQ
LGLKDFSPKSQDAVALQQIK
ERGALPMIDRGDIRQAIDRC
SNIWASLPGAGYGQFEHKAD
SLIAKFKEAGGTVREIDV
