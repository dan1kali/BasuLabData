# Requirements
- MATLAB 2024b or later (for plotting) + Signal Processing Toolbox
- Chronux functions, generateCrossValInd.m, permutationTest.m added to path

# Usage

### 1. <u>Place all patient data in folder</u>
**a/b) Time series EEG**
- Fieldtrip notation
- use folder called "ab_inputEEGData"
- each patient should have their own .mat file with data saved in ft_data3_filt_rs
- and other variables

**c) Preprocessed Time series spectral power**
- Baselined, channels selected for conflict modulation
- use folder called "c_inputPowerData"
- each patient should have their own FOLDER with data saved in: normPow_HighGamma, normPow_Theta etc
- and other variables

### 2. <u>Generate Intermediate Data</u>
**a_chronux.m - Generate power data, features, and other necessary indices with chronux toolbox**
- Preprocess EEG, generate powerdata with Chronux + other files (non z-scored features, time, C/I trials, etc)  with preprocess function
- Generate conflict modulated channels + z scored power data with conflictModAnalysis function
- All files will save in auto-generated outputData folders for each patient

**b_wavelet.m - Generate power using pre-processed EEG data, extract features with wavelet**
- Manual import of preprocessed EEG + artifact trials
- Generate power data using wavelet with preprocess function
- Implements baselining, conflict modulation analysis, feature extraction with conflictModAnalysis function
- (Wavelet extraction, updated simplified baselining, and eliminated responsive channel analysis)
- saves in folder b_outputPowerData

**c_features.m - Extract features only**
- Manual import of normalized power data
- Feature extraction only
- saves in folder c_outputPowerData

### 3. <u>Apply decoding and plot</u>
**decode.m**
- Update patient list and group conditions
- Load output data saved in patient save folders obtained in step 2
- Use nBar = 1 to create one bar grouping for aggregated multi-patient data seperated by condition
- Or nBar = length(patients) for bar groupings seperated by patient AND condition

# Obtain Plots Only

**Use plots.m**
- generates plots for accuracies by region, grouped region, and patient