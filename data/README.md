## Data Description

### mgtstMetadata
- `seq_lab`: laboratory where the sequencing library was generated and sequencing run was performed, JHU or NIST. 	 
- `seq_run`: sequencing run number, 1 or 2, sequencing runs at the two sequencing labs are independent.  	
- `sample_id`: combined `seq_lab`, `seq_run`, `pcr_16S_plate`, and `pos` variables.  
- `biosample_id`: vaccine trail participant id.    	
- `t_fctr`: sample titration number 1-5, 10, and 15 are titrations, 0 is the unmixed post-exposure sample, and 20 is the unmixed pre-exposure sample.   
- `pcr_16S_plate`: 16S PCR plate replicate number, 1 or 2
- `pos`: position in 96 well plate
