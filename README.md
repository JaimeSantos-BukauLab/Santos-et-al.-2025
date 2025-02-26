Python scripts used for binding period detection.
Additional data analysis of the ribosome profiling datasets (i.e: Single gene and metagene analysis) was performed using the RiboSeqTools package [available at https://github.com/ilia-kats/RiboSeqTools]


-1_Convert_h5file_to_pickle.py uses h5files as an input to generate pickle dictionaries for the following steps.

-2_Confidence_interval_calc.py calculates confidence intervals based on the Selectome (IP) and Total translatomes (TT) using defined windows

-3_Get_binding_periods.py Extract binding periods using the heuristic rule described in the method section of the article

-00_Sequence_human_proteins.pick and 00_Localization_human_proteins.pick are support dictionaries containing human protein sequences and protein localization annotation.
