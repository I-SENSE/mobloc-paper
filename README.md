# Mobintel: CSI-based Localization

This repository contains the source code for MobLoc - a CSI-based fingerprinting method. This method has been described in a paper, "MobLoc: CSI-based Location Fingerprinting with MUSIC". The paper has been accepted by [IEEE Journal of Indoor and Seamless Positioning and Navigation (J-ISPIN)](https://ieee-jispin.org/).

Authors:
* [Stephan Mazokha](https://scholar.google.com/citations?user=tCc9M3EAAAAJ&hl=en&oi=ao)
* [Fanchen Bao](https://scholar.google.com/citations?user=LAAw4LMAAAAJ&hl=en)
* [George Sklivanitis](https://scholar.google.com/citations?user=kNTCxpEAAAAJ&hl=en)
* [Jason O. Hallstrom](https://www.fau.edu/engineering/directory/faculty/hallstrom/)

For further details, reach out to [I-SENSE at Florida Atlantic University](https://isense.fau.edu/).

## How to run this code

0. Install Matlab;
1. Download DLoc dataset [from this website](https://wcsng.ucsd.edu/wild);
2. Save the channel MAT files in `./data/channels`, split indexes in `./data/split_idx`;
3. To reproduce the CDF chart from Figure 9, launch `./mobloc/mobloc_simple_env.m`;
4. To reproduce the CDF chart from Figure 10, launch `./mobloc/mobloc_complex_env.m`.

## Additional resources

* DLoc repository: [GitHub](https://github.com/ucsdwcsng/DLoc_pt_code)
* Mobintel website (core I-SENSE project): [Website](https://www.mobintel.org)
