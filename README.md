**testRun2.m - apply decoding and plot**
- plotting requires some manuvering and uncommenting 

**extractChronuxFeatures4.m - most recent**
- updates for new data + parcellations + preprocessing
- will update for file format to move away from structure + paralellizing

- must run preprocessing step first to generate features file, then conflict modulation step

old:

extractChronuxFeatures.m - original with band power not z score

extractChronuxFeatures2.m - with z score

extractChronuxFeatures3.m - changing the way the band power is cut
- Add preprocessing support for their data

extractChronuxFeatures3_1.m - small updates
- Extract power for all channels at once + other changes
