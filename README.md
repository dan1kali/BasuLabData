# Requirements
- MATLAB 2024b or later (for plotting) + Signal Processing Toolbox
- Chronux functions, generateCrossValInd.m, permutationTest.m added to path

# Usage

### 1) Place all patient data in folder called "patientData"

### 2) Use extractChronuxFeatures4.m
**generate power data, features, and all other nexessary files**
- Use preprocessing function to generate powerdata + other files (non z-scored features, time, C/I trials, etc)
- Use conflictModAnalysis function to generate conflict modulated channels + z scored power data
- All files will save in auto-generated outputData folders for each patient

### 3) testRun3.m
**apply decoding and plot**
- update patient list and group conditions
- use nBar = 1 to create one bar grouping for aggregated multi-patient data seperated by condition
- or nBar = length(patients) for bar groupings seperated by patient AND condition
- Will load from data saved in patient save folders using previous file



## Other Files

**extractChronuxFeatures5.m**
- for use with wavelet extraction, updated simplified baselining, and eliminated responsive channel analysis

**extractChronuxFeatures6.m**
- manual import of power data, manual import of conflict modulated channels