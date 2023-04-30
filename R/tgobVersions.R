#' @title Subsample versions (TO, TS and ssosX_\emph{species}) from TGOB
#'
#' @param p sf/sfc (POINT/MULTIPOINT): point occurrences (TGOB)
#' @param r RasterLayer: Template raster. If provided, all stats are calculated per raster pixel (cell/square).
#' @param species character: column name with species name
#' @param TS.n integer or numeric: int or ratio (< 1, \code{cumsum} \emph{p}) to select for TOP species
#' @param observers character: column name with observers' names (expected: "John Doe" / "Doe J." / "John Doe, Jane Maria Moe" / "Doe J., Moe J. M." / ...)
#' @param TO.n integer or numeric: int or ratio (< 1, \code{cumsum} \emph{p}) to select for TOP observers
#' @param observersRemoveSingleName logical: remove single word (name) observers
#' @param observersRemoveSingleOccurrence integer: remove observers with single (or specified number) species occurrence from ssos. Can be too restrictive for very rare species, where observers have small chance to re-observe. On the other hand, can be useful to restrict very active observers that record some species only purposefully.
#' @param quantileObserversThreshold numeric: (0-1) (centile 0.01, decile 0.10, ...) threshold to remove individual lower outlier observers with "suspicious" (purposefull) pattern, calculated as \emph{observational ratio} (Σ\emph{focus species} / Σ\emph{rest of species}) compared to overall observers \emph{observational ratio}. In other words: comparing (per species, \emph{ssos}?) each observers' profile (sums of species occurerrences) to overall observers' profile and then remove outliers.
#' @param TSAO.n integer or numeric: int or quantile (< 1, upper quartile 0.75, ...) to select for TOP species by Area Overlap (overlaping occupied cells with focus species). Can be done only with inserted \emph{r}!
#' @param TSAO.min numeric (0-1): minimum ratio overlap to threshold before \emph{TSAO.n} is applied
#' @param badWordsSpecies character vector: if species name containing unwanted strings, remove such rows
#' @param badWordsObservers character vector: if observers name containing unwanted strings, remove such rows
#' @param crs integer: Force crs.
#' @param prefix character: prefix for newly added columns to \emph{p} with binary (0/1) sign that links affiliation to version (to select it easily later)
#'
#' @return named list: \emph{t}: sf (POINT/MULTIPOINT): inserted \emph{p} with new (0/1) columns: TO, TS, ssosX_\emph{species}; \emph{report}: TO+TS stats used to threshold selected top O+S
#'
#' @export

tgobVersions <- function(p, r = NA, species = "species", TS.n = 0.2, observers = NA, TO.n = 0.2, observersRemoveSingleName = TRUE, observersRemoveSingleOccurrence = 0, quantileObserversThreshold = 0.01, TSAO.n = 0.8, TSAO.min = 0.5, badWordsSpecies = NA, badWordsObservers = NA, crs = NA, prefix = "nc_") {
    badWords <- "_badWords"

    out <- list("t" = NA, "report" = NA)

    if (!is(p, "sf") & !(sf::st_geometry_type(p, by_geometry = FALSE) == "POINT")) {
        stop("p: only sf/sfc (POINT/MULTIPOINT) allowed!")
        return(out)
    }

    if (is.integer(crs)) {
        sf::st_crs(p) <- paste0("epsg:", crs)
        if (is(r, "RasterLayer")) sf::st_crs(r) <- paste0("epsg:", crs)
    } else {
        if (is(r, "RasterLayer") & !(sf::st_crs(p)$proj4string == sf::st_crs(r)$proj4string)) {
            stop("p and r crs are not equal!")
            return(out)
        }
        crs <- raster::crs(r)
    }

    p %<>% ungroup() %>% mutate(uid = row_number())

    # # # # # # # # # #
    # prepare ssos  + prefilter badWords
    # # # # # # # # # #

    if (!is.na(observers)) {
        # save original observers' names
        p %<>% mutate(!!paste0(observers, "_orig") := !!sym(observers))

        # remove non-latin characters + lowercase
        p %<>% mutate(!!observers := tolower(stringi::stri_trans_general(str = !!sym(observers), id = "Latin-ASCII")))

        # make two words from dot separated strings
        p %<>% mutate(!!observers := str_replace(!!sym(observers), "\\.", "\\. "))
        # trim spaces
        p %<>% mutate(!!observers := trimws(!!sym(observers)))

        # mark single word observers
        p %<>% mutate(!!paste0(observers, "_single") := ifelse(stringi::stri_count_fixed(!!sym(observers), " ") == 0, 1, 0))

        # remove single word observers
        if (observersRemoveSingleName) {
            p %<>% filter(!!sym(paste0(observers, "_single")) == 0)
        }

        if (!is.na(badWordsObservers)) {
            # mark bad words observers
            p %<>% mutate(!!paste0(observers, badWords) := ifelse(str_detect(!!sym(observers), paste(badWordsObservers, collapse = "|")), 1, 0))
        } else {
            p %<>% mutate(!!paste0(observers, badWords) := 0)
        }
        p %<>% filter(!!sym(paste0(observers, badWords)) == 0)
    }

    if (!is.na(badWordsSpecies)) {
        # mark bad words species
        p %<>% mutate(!!paste0(species, badWords) := ifelse(str_detect(!!sym(species), paste(badWordsSpecies, collapse = "|")), 1, 0))
    } else {
        p %<>% mutate(!!paste0(species, badWords) := 0)
    }

    p %<>% filter(!!sym(paste0(species, badWords)) == 0)

    if (is(r, "RasterLayer")) {
        #
        # per species, (observer) and pixel
        #

        v.temp <- raster::extract(r, sf::st_coordinates(p), cellnumbers = TRUE)
        p %<>% tibble::add_column(cellNumber = v.temp[, 1])
        p %<>% filter(!is.na(cellNumber))

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(!!sym(observers), cellNumber, uid, -geometry)) %>%
                ungroup() %>%
                group_by(!!sym(observers), cellNumber) %>%
                slice_head(n = 1)
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(!!sym(species), cellNumber, uid, -geometry)) %>%
            ungroup() %>%
            group_by(!!sym(species), cellNumber) %>%
            slice_head(n = 1)
    } else {
        #
        # per species and (observer)
        #

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(!!sym(observers), uid, -geometry)) %>%
                group_by(!!sym(observers))
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(!!sym(species), uid, -geometry)) %>%
            group_by(!!sym(species))
    }

    # # # # # # # # # #
    # TO
    # # # # # # # # # #

    if (!is.na(observers)) {
        # per observer stats

        uid.total <- nrow(p.temp.o)
        p.TO.stat <- p.temp.o %>%
            ungroup() %>%
            group_by(!!sym(observers)) %>%
            summarise(
                uid.n = n_distinct(uid),
                uid.ratio = n_distinct(uid) / uid.total
            ) %>%
            arrange(desc(uid.n))

        if (TO.n >= 1) {
            # exact number of TOP observers
            TO.total <- p.temp.o %>% n_distinct(as.character(sym(observers)))
            if (TO.n > TO.total) {
                stop(paste0("Can't get more TOP observers than total number: ", TO.total))
            } else {
                p.TO.unique <- p.TO.stat %>%
                    arrange(desc(uid.n)) %>%
                    slice_head(n = TO.n)
            }
        } else {
            p.TO.stat %<>% add_column(cumsum = cumsum(p.TO.stat$uid.ratio))
            p.TO.unique <- p.TO.stat %>% filter(cumsum <= TO.n)
        }

        p.TO.unique <- unique(as.vector(unlist(p.TO.unique %>% dplyr::select(!!sym(observers)))))
        # mark top X observers
        p %<>% mutate(!!paste0(prefix, "TO") := ifelse(!!sym(observers) %in% p.TO.unique, 1, 0))
    }

    # # # # # # # # # #
    # TS
    # # # # # # # # # #

    # per species stats

    uid.total <- nrow(p.temp.s)
    p.TS.stat <- p.temp.s %>%
        ungroup() %>%
        group_by(!!sym(species)) %>%
        summarise(
            uid.n = n_distinct(uid),
            uid.ratio = n_distinct(uid) / uid.total
        ) %>%
        arrange(desc(uid.n))

    if (TS.n >= 1) {
        # exact number of TOP species
        TS.total <- p.temp.s %>% n_distinct(as.character(sym(species)))
        if (TS.n > TS.total) {
            stop(paste0("Can't get more TOP species than total number: ", TS.total))
        } else {
            p.TS.unique <- p.TS.stat %>%
                arrange(desc(uid.n)) %>%
                slice_head(n = TS.n)
        }
    } else {
        p.TS.stat %<>% add_column(cumsum = cumsum(p.TS.stat$uid.ratio))
        p.TS.unique <- p.TS.stat %>% filter(cumsum <= TS.n)
    }

    p.TS.unique <- unique(as.vector(unlist(p.TS.unique %>% dplyr::select(!!sym(species)))))
    # mark top X species
    p %<>% mutate(!!paste0(prefix, "TS") := ifelse(!!sym(species) %in% p.TS.unique, 1, 0))

    # # # # # # # # # #
    # ssos - add species/version columns (0/1)
    # # # # # # # # # #


    species.unique <- unique(as.vector(unlist(as_tibble(p) %>% dplyr::select(!!sym(species)))))

    # ssos versions
    p.t <- as_tibble(p) %>% dplyr::select(-geometry)
    TSAO.min.species <- c()
    TSAO.top <- list()

    # species
    for (sp in species.unique) {
        ssos.temp <- NA
        message(sp)
        p.t.sp <- p.t %>% filter(!!sym(species) == sp)

        if (nrow(p.t.sp) < 1) {
            # zero species occurrences in ssos
            next
        }

        # select observers with at least one occ of selected species
        observers.unique <- unique(as.vector(unlist(p.t.sp %>% dplyr::select(!!sym(observers)))))
        ssos.temp <- p.t %>% filter(!!sym(observers) %in% observers.unique)

        # count occs per species and observers
        ssos.temp.stat <- ssos.temp %>%
            ungroup() %>%
            group_by(!!sym(species), !!sym(observers)) %>%
            summarise(uid.n = n_distinct(uid))

        # sum occs per observers
        ssos.temp.stat.observers <- ssos.temp.stat %>%
            ungroup() %>%
            group_by(!!sym(observers)) %>%
            summarise(observers.n = sum(uid.n))

        # sum occs of selected species per observers
        ssos.temp.stat.species <- ssos.temp.stat %>%
            ungroup() %>%
            filter(!!sym(species) == sp) %>%
            group_by(!!sym(observers)) %>%
            summarise(species.n = sum(uid.n))

        # ratio (selected species occs per total occs) per observers
        ssos.temp.ratio <- ssos.temp.stat.observers %>%
            left_join(ssos.temp.stat.species, by = as.character(sym(observers))) %>%
            mutate(ratio = species.n / observers.n)

        observersRemoveSingleOccurrence <- as.numeric(observersRemoveSingleOccurrence)
        if (observersRemoveSingleOccurrence > 0) {
            # remove single
            ssos.temp.ratio %<>% filter(species.n > observersRemoveSingleOccurrence)
            observers.unique <- unique(as.vector(unlist(ssos.temp.ratio %>% dplyr::select(!!sym(observers)))))
        }

        quantileObserversThreshold <- as.numeric(quantileObserversThreshold)
        if (quantileObserversThreshold > 0 & quantileObserversThreshold < 1) {
            # remove quantile (and/or singles before)
            # remove "outlier" observers with suspicious ratio (lower than input quantile)
            ssos.temp.ratio.quantile <- unname(stats::quantile(ssos.temp.ratio$ratio, probs = quantileObserversThreshold, na.rm = TRUE))
            observers.unique <- unique(as.vector(unlist(ssos.temp.ratio %>% filter(ratio > ssos.temp.ratio.quantile[1]) %>% dplyr::select(!!sym(observers)))))
        }

        p %<>% mutate(!!paste0(prefix, "ssos", "_", sp) := ifelse(!!sym(observers) %in% observers.unique, 1, 0))

        # # # # # # # # # #
        # TSAO
        # # # # # # # # # #
        if (is(r, "RasterLayer")) {
            ss.total <- length(unique(p.t.sp$cellNumber))

            p.TSAO.stat <- p.temp.s %>%
                ungroup() %>%
                group_by(!!sym(species)) %>%
                filter(cellNumber %in% unique(p.t.sp$cellNumber)) %>%
                summarise(
                    cells.shared = n_distinct(cellNumber),
                    cells.shared.ratio = n_distinct(cellNumber) / ss.total
                ) %>%
                arrange(desc(cells.shared))

            # remove focus species (100 % overlap...)
            # p.TSAO.stat %<>% filter(!!sym(species) != sp)

            if (nrow(p.TSAO.stat) > 0) {
                if (TSAO.n >= 1) {
                    # exact number of TOP species
                    p.TSAO.unique <- p.TSAO.stat %>% slice_head(n = TSAO.n)
                } else {
                    if (nrow(p.TSAO.stat %>% filter(cells.shared.ratio > TSAO.min)) >= 4) {
                        # left only >50 % overlaping species with more than 3 species
                        p.TSAO.stat.temp <- p.TSAO.stat %>% filter(cells.shared.ratio > TSAO.min)
                        tsao.quantile <- NA
                        tsao.quantile <- unname(stats::quantile(p.TSAO.stat.temp$cells.shared, probs = TSAO.n, na.rm = TRUE))

                        if (nrow(p.TSAO.stat.temp %>% filter(cells.shared > tsao.quantile)) >= 4) {
                            p.TSAO.unique <- p.TSAO.stat.temp %>% filter(cells.shared > tsao.quantile)
                        } else {
                            message(paste0("Less than 3 ovelaping species with >", TSAO.min, " overlap and quantile ", TSAO.n, ". Left 3 top species"))
                            TSAO.min.species <- c(TSAO.min.species, sp)
                            p.TSAO.unique <- p.TSAO.stat %>% slice_head(n = 4)
                        }
                    } else {
                        message(paste0("Less than 3 ovelaping species with >", TSAO.min, " overlap. Left 3 top species"))
                        TSAO.min.species <- c(TSAO.min.species, sp)
                        p.TSAO.unique <- p.TSAO.stat %>% slice_head(n = 4)
                    }
                    TSAO.top[[sp]] <- p.TSAO.unique
                }
            } else {
                message("No overlaping species...")
                next
            }

            p.TSAO.unique <- unique(as.vector(unlist(p.TSAO.unique %>% dplyr::select(!!sym(species)))))

            # mark top X  AO species
            p %<>% mutate(!!paste0(prefix, "TSAO", "_", sp) := ifelse(!!sym(species) %in% p.TSAO.unique, 1, 0))
        }
    }

    return(list(
        "t" = p,
        "report" = list(
            "TS" = p.TS.stat, "TS.n" = TS.n,
            "TO" = p.TO.stat, "TO.n" = TO.n,
            "TSAO" = TSAO.top, "TSAO.n" = TSAO.n, "TSAO.min.species" = TSAO.min.species
        )
    ))
}
