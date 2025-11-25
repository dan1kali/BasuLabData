# Requirements
- MATLAB 2024b or later (for plotting) + Signal Processing Toolbox
- Chronux functions, generateCrossValInd.m, permutationTest.m added to path

# Usage

### 1. <u>Place all patient data in folder</u>
**a) Raw time series EEG**
- Fieldtrip notation
- use folder called "patientData"

**b) Time series spectral power**
- Unprocessed
- use folder called "patientData"

**c) Preprocessed Time series spectral power**
- Baselined, channels selected for conflict modulation
- use folder called "patientPowerData"

### 2. <u>Generate Intermediate Data</u>
**a) Generate power data, features, and all other necessary indices: Use extractChronuxFeatures4.m**
- Use preprocessing function to generate powerdata + other files (non z-scored features, time, C/I trials, etc)
- Use conflictModAnalysis function to generate conflict modulated channels + z scored power data
- All files will save in auto-generated outputData folders for each patient

**b) Process pre-existing power data and extract features - extractWaveletFeatures.m**
- Manual import of preprocessed power data + artifact trials
- Implements baselining, conflict modulation analysis, feature extraction
- (Wavelet extraction, updated simplified baselining, and eliminated responsive channel analysis)

**c) Extract features only - extractFeaturesOnly.m**
- Manual import of already normalized power data
- Feature extraction only

### 3. <u>Apply decoding and plot</u>
**Use testRun4.m**
- Update patient list and group conditions
- Use nBar = 1 to create one bar grouping for aggregated multi-patient data seperated by condition
- Or nBar = length(patients) for bar groupings seperated by patient AND condition
- Will load from data saved in patient save folders using previous file

# Obtain Plots Only

**Use plots.m**
- generates plots for accuracies by region, grouped region, and patient