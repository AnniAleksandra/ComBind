Use R software version > 3.5.2

*NOTE: Training ComBind is slow and requires gigabytes of RAM with high-throughput sequencing data
*Required R-packages: ShortRead, randomForest, reshape2, pROC

##### 1. Place training set DNA
Place positive training sequences to: DNA/CAP-SELEX-Positive-Training/
Place negative training sequences to: DNA/CAP-SELEX-Negative-Training/

1. Name the data sets in training folders with the same name
2. Ideally, positive and negative data sets in training should be balanced in size
3. Format: fasta.gz
4. Sequence length is at least 25 nucleic acids, preferably 40 nucleic acids, while upper limit is not defined
5. All sequences should have the same length
6. Sequences with ambigious letters/nucleotides are removed from random forest training

##### 2. Train ComBind
1. Run from the command line in ComBind folder:
Rscript ./CODE/Train_ComBind.R DATASET_NAMES_IN_DNA_FOLDERS TF1 TF2

*NOTE: Make sure that 'PWM_HT-SELEX' folder include PWMs for TF1 and TF2 (names accordingly)

2. Example for Alx4-Tbx21 using a subset of CAP-SELEX DNA sequences:
Rscript ./CODE/Train_ComBind.R ALX4_TBX21_SUBSET ALX4 TBX21

##### 3. Place test set DNA
Place sequences you wish to test with ComBind to: DNA/SCORE-DNA/

1. Format: fasta.gz
2. Sequence length is at least 25 nucleic acids, preferably 40 nucleic acids, while upper limit is not defined (sequence length can be different from training sequence length)
3. All testing sequences must have the same length
4. Testing sequences should not entail ambigious letters (only "A", "T", "C" and "G" are allowed)

##### 4. Test ComBind
1. Run from the command line in ComBind folder:
Rscript ./CODE/Test_ComBind.R TESTSET_NAME_IN_DNA_FOLDER RF_NAMES_IN_RF_FOLDER_WITHOUT_RF_NUMBER_PREFIXES

2. Example for Alx4-Tbx21 using a subset of CAP-SELEX DNA sequences:
Rscript ./CODE/Test_ComBind.R ALX4_TBX21 ALX4_TBX21_SUBSET

