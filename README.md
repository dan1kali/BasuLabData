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
**a) Generate power data, features, and other necessary indices: a_chronux.m**
- Preprocess EEG, generate powerdata with Chronux + other files (non z-scored features, time, C/I trials, etc)  with preprocess function
- Generate conflict modulated channels + z scored power data with conflictModAnalysis function
- All files will save in auto-generated outputData folders for each patient

**b) Generate power using wavelet using pre-processed data, extract features: b_wavelet.m**
- Manual import of preprocessed EEG + artifact trials
- Generate power data using wavelet with preprocess function
- Implements baselining, conflict modulation analysis, feature extraction with conflictModAnalysis function
- (Wavelet extraction, updated simplified baselining, and eliminated responsive channel analysis)

**c) Extract features only: c_features.m**
- Manual import of normalized power data
- Feature extraction only

### 3. <u>Apply decoding and plot</u>
**decode.m**
- Update patient list and group conditions
- Use nBar = 1 to create one bar grouping for aggregated multi-patient data seperated by condition
- Or nBar = length(patients) for bar groupings seperated by patient AND condition
- Will load from data saved in patient save folders using previous file

# Obtain Plots Only

**Use plots.m**
- generates plots for accuracies by region, grouped region, and patient