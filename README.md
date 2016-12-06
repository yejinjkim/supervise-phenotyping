# Discriminative and Distinct Phenotyping by Nonnegative Tensor Factorization 
DEC06 2016 Ver0.0

This is a repository of code for Discriminative and Distinct Phenotyping by Nonnegative Tensor Factorization. 

## REQUIREMENTS
This code requires the following:
- Matlab (www.mathworks.com)
- Tensor Toolbox Version 2.6 (http://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html)


## DATA FORMAT
- count.csv: co-occurrence count of diagnosis and prescription for each subject. 

		count, subject_id, diagnosis_id, prescription_id
		1, 35, 16, 41
		2, 35, 16, 68
		...

- label.csv : subject id of clinical outcome (eg. mortality, readmission) is 1 (eg. deceased, readmitted)

		subject_id
		2
		5
		10
		45
		...
- similarities.csv: pairwise similarity matrix of diagnoses and prescriptions.

							d_diabetes, d_hypertension, ..., p_propofol, p_morphine_sulfate, ...
		d_diabetes				1,				0.3,				
		d_hypertension			0.3, 			1,					0.2
		...
		p_propofol
		p_morphine_sulfate
		...
## RUN
1. Prepare CSV files: count.csv, label.csv, similarities.csv
2. Run main.m

