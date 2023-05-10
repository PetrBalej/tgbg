# *tgbg*: Target-Group Background Generator

*Target-Group (Occurrences) Background (and Bias Raster) Generator*

## Variants:
1) **TGOB**
2) **TO**: Top X Observers subsample of TGOB
3) **TS**: Top X Species subsample of TGOB
4) **TSAO**\*: Top X Species by Area Overlap subsample of TGOB
5) **sso**\*: *Shared Species-Observers* subsample of previous versions: TGOBsso, TOsso, TSsso, TSAOsso

\* In order to standardise the intensity of sampling effort with respect to individual species, we propose two new concepts: 

ad 4) Focus species specific modification of 3). Select for TOP X species by Area Overlap (most overlaping occupied cells with focus species). We assume that it link (more observed) spatially related co-observed species with appropriate sampling effort intensity

ad 5) "*Shared Species-Observers* subsample of previous versions". The principle is to subsample from previous variants (1-4) only those observers who have already observed the selected species at least once (or more). This simple rule can be implemented at the query level of the species occurrence database by matching the names of observers and filtering out observers with zero records of selected species. We assume that it removes potentially detrimental noise from sampling effort from observers who are not sure that they are able to determine the given species or do not map it purposefully. 

![Variants: TGOB, TO, TS, TSAO, TGOBsso, TOsso, TSsso, TSAOsso](/files/diagram.png)

## Sources (basics):

Dudík, M., Phillips, S. and Schapire, R.E., 2005. Correcting sample selection bias in maximum entropy density estimation. Advances in neural information processing systems, 18. <https://dl.acm.org/doi/10.5555/2976248.2976289>

Phillips, S.J., Dudík, M., Elith, J., Graham, C.H., Lehmann, A., Leathwick, J. and Ferrier, S., 2009. Sample selection bias and presence‐only distribution models: implications for background and pseudo‐absence data. Ecological applications, 19(1), pp.181-197. <https://doi.org/10.1890/07-2153.1>

## Sources (advanced):

Ranc, N., Santini, L., Rondinini, C., Boitani, L., Poitevin, F., Angerbjörn, A. and Maiorano, L., 2017. Performance tradeoffs in target‐group bias correction for species distribution models. Ecography, 40(9), pp.1076-1087. <https://doi.org/10.1111/ecog.02414>

Vollering, J., Halvorsen, R., Auestad, I. and Rydgren, K., 2019. Bunching up the background betters bias in species distribution models. Ecography, 42(10), pp.1717-1727. <https://doi.org/10.1111/ecog.04503>

Botella, C., Joly, A., Monestiez, P., Bonnet, P. and Munoz, F., 2020. Bias in presence-only niche models related to sampling effort and species niches: Lessons for background point selection. PLoS One, 15(5), p.e0232078. <https://doi.org/10.1371/journal.pone.0232078>

Komori, O., Eguchi, S., Saigusa, Y., Kusumoto, B. and Kubota, Y., 2020. Sampling bias correction in species distribution models by quasi-linear Poisson point process. Ecological Informatics, 55, p.101015. <https://doi.org/10.1016/j.ecoinf.2019.101015>

Chauvier, Y., Zimmermann, N.E., Poggiato, G., Bystrova, D., Brun, P. and Thuiller, W., 2021. Novel methods to correct for observer and sampling bias in presence‐only species distribution models. Global Ecology and Biogeography, 30(11), pp.2312-2325. <https://doi.org/10.1111/geb.13383>