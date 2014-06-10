CAPPED Dipeptides ACE-X-NME
===========================

The directory population_data contains three files containing estimated populations of the alpha, beta, and PPII 
conformational states.  The datasets are from the Baldwin and the Sosnick papers.

baldwin_table1.csv contains the infrared column from table 1.  When multiple pH values were available, we
chose the pH 4.9 values.  

baldwin_table1.csv is table 1 from the 2006 baldwin paper, containing scalar couplings of the capped dipeptides at 4.9 pH

Because HIS has a pKA near 5.0, our code manually HARD-CODES the value for 2.9 pH.  This allows us to use a fixed protonation state
rather than a constant pH simulation.  
