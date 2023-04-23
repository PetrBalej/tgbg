# tgbg (Target-Group Background Generator)

*Target-Group (Occurrences) Background (and Bias Raster) Generator*

## Variants:
1) **TGOB**
2) **TO**: Top X Observers
3) **TS**: Top X Species
4) **ssos**\*: *Shared Species-Observers Subsample* of previous versions (ssosTGOB, ssosTO, ssosTS)

\* In order to standardise the intensity of sampling effort with respect to individual species, we propose a new concept: "*Shared Species-Observers Subsample*" (**ssos**). The principle is to subsample from TGOB only those observers who have already observed the selected species at least once (or more). This simple rule can be implemented at the query level of the species occurrence database by matching the names of observers and filtering out observers with zero records of selected species. We assume that it removes potentially detrimental noise from sampling effort from observers who are not sure that they are able to determine the given species or do not map it purposefully. 
