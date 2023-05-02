# tgbg (Target-Group Background Generator)

*Target-Group (Occurrences) Background (and Bias Raster) Generator*

## Variants:
1) **TGOB**
2) **TO**: Top X Observers subsample of TGOB
3) **TS**: Top X Species subsample of TGOB
4) **TSAO**\*: Top X Species by Area Overlap subsample of TGOB
5) **sso**\*: *Shared Species-Observers* subsample of previous versions: ssoTGOB, ssoTO, ssoTS, ssoTSAO

\* In order to standardise the intensity of sampling effort with respect to individual species, we propose two new concepts: 

ad 4) Focus species specific modification of 3). Select for TOP X species by Area Overlap (most overlaping occupied cells with focus species). We assume that it link (more observed) spatially related co-observed species with appropriate sampling effort intensity

ad 5) "*Shared Species-Observers* subsample of previous versions". The principle is to subsample from previous variants (1-4) only those observers who have already observed the selected species at least once (or more). This simple rule can be implemented at the query level of the species occurrence database by matching the names of observers and filtering out observers with zero records of selected species. We assume that it removes potentially detrimental noise from sampling effort from observers who are not sure that they are able to determine the given species or do not map it purposefully. 
