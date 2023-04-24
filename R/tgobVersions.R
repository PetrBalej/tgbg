#' @title Subsample versions (TO, TS, ssosTO, ssosTS, ssosTGOB) from TGOB
#'
#' @param p sf/sfc (POINT/MULTIPOINT): Points
#' @param r RasterLayer: Template raster. If provided, all stats are calculated per raster pixel (cell/square).
#' @param species character: column name with species name
#' @param TS.n integer or numeric: Ratio (< 1) or int to select for TOP species
#' @param observers character: column name with observers' names (expected: "John Doe" / "Doe J." / "John Doe, Jane Maria Moe" / "Doe J., Moe J. M." / ...)
#' @param TO.n integer or numeric: Ratio (< 1) or int to select for TOP observers
#' @param observersRemoveSingle logical: remove single word observers
#' @param badWordsSpecies character vector: observers containing unwanted strings
#' @param badWordsObservers character vector: observers containing unwanted strings
#' @param crs integer: Force crs.
#'
#' @return sf (POINT/MULTIPOINT): inserted \emph{p} with new (0/1) columns: TO, TS, ssosTGOB_\emph{species}, ssosTO_\emph{species}, ssosTS_\emph{species}

tgobVersions <- function(p, r = NA, species = "species", TS.n = 0.1, observers = NA, TO.n = 0.1, observersRemoveSingle = TRUE, badWordsSpecies = NA, badWordsObservers = NA, crs = NA) {
    badWords <- "_badWords"
    prefix <- "nc_" # new columns with versions to detect

    if (!is.na(r) & !is(r, "RasterLayer")) {
        stop("r: only RasterLayer allowed!")
        return(NA)
    }

    if (is(p, "sf") && sf::st_geometry_type(rp, by_geometry = FALSE) == "POINT") {
        stop("p: only sf/sfc (POINT/MULTIPOINT) allowed!")
        return(NA)
    }

    if (is.integer(crs)) {
        sf::st_crs(p) <- paste0("epsg:", crs)
        sf::st_crs(r) <- paste0("epsg:", crs)
    } else {
        if (sf::st_crs(p)$proj4string == sf::st_crs(r)$proj4string) {
            stop("p and r crs are not equal!")
            return(NA)
        }
        crs <- raster::crs(r)
    }

    p %<>% ungroup() %>% mutate(uid = row_number())

    # # # # # # # # # #
    # prepare ssos  + prefilter badWords
    # # # # # # # # # #

    if (!is.na(observers)) {
        # save original observers' names
        p %<>% mutate(sym(paste0(observers, "_orig")) := sym(observers))
        # remove non-latin characters + lowercase
        p %<>% mutate(sym(observers) := tolower(stringi::stri_trans_general(str = sym(observers), id = "Latin-ASCII")))
        # make two words from dot separated strings
        p %<>% mutate(sym(observers) := str_replace(sym(observers), "\\.", "\\. "))
        # trim spaces
        p %<>% mutate(sym(observers) := trimws(sym(observers)))

        # mark single word observers
        p %<>% mutate(sym(paste0(observers, "_single")) := ifelse(stringi::stri_count_fixed(sym(observers), " ") == 0, 1, 0))

        # remove single word observers
        if (observersRemoveSingle) {
            p %<>% filter(sym(paste0(observers, "_single")) == 0)
        }

        if (!is.na(badWordsObservers)) {
            # mark bad words observers
            p %<>% mutate(sym(paste0(observers, badWords)) := ifelse(str_detect(sym(observers), paste(badWordsObservers, collapse = "|")), 1, 0))
        } else {
            p %<>% mutate(sym(paste0(observers, badWords)) := 0)
        }
        p %<>% filter(sym(paste0(observers, badWords)) == 0)
    }

    if (!is.na(badWordsSpecies)) {
        # mark bad words species
        p %<>% mutate(sym(paste0(species, badWords)) := ifelse(str_detect(sym(species), paste(badWordsSpecies, collapse = "|")), 1, 0))
    } else {
        p %<>% mutate(sym(paste0(species, badWords)) := 0)
    }

    p %<>% filter(sym(paste0(species, badWords)) == 0)

    if (!is.na(r)) {
        #
        # per species, (observer) and pixel
        #

        v.temp <- raster::extract(r, sf::st_coordinates(p))
        p$cellNumber <- v.temp
        p %<>% filter(!is.na(cellNumber))

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(sym(observers), cellNumber, uid, -geometry)) %>%
                ungroup() %>%
                group_by(sym(observers), cellNumber) %>%
                slice_head(n = 1)
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(sym(species), cellNumber, uid, -geometry)) %>%
            ungroup() %>%
            group_by(sym(species), cellNumber) %>%
            slice_head(n = 1)
    } else {
        #
        # per species and (observer)
        #

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(sym(observers), uid, -geometry)) %>%
                group_by(sym(observers))
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(sym(species), uid, -geometry)) %>%
            group_by(sym(species))
    }

    # # # # # # # # # #
    # TO
    # # # # # # # # # #

    if (!is.na(observers)) {
        # per observer stats

        uid.total <- nrow(p.temp.o)
        p.TO.stat <- p.temp.o %>%
            summarise(
                uid.n = n_distinct(uid),
                uid.ratio = n_distinct(uid) / uid.total
            ) %>%
            arrange(desc(uid.n))

        if (TO.n >= 1) {
            # exact number of TOP observers
            TO.total <- p.temp.o %>% n_distinct(sym(observers))
            if (TO.n > TO.total) {
                stop(paste0("Can't get more TOP observers than total number: ", TO.total))
            } else {
                p.TO.unique <- p.TO.stat %>%
                    arrange(desc(uid.n)) %>%
                    slice_head(n = TO.n)
            }
        } else {
            p.TO.unique <- p.TO.stat %>% filter(uid.ratio <= TO.n)
        }

        p.TO.unique <- unique(as.vector(unlist(p.TO.unique)))
        # mark top X observers
        p %<>% mutate(sym(paste0(prefix, "TO")) := ifelse(sym(observers) %in% p.TO.unique, 1, 0))
    }

    # # # # # # # # # #
    # TS
    # # # # # # # # # #

    # per species stats

    uid.total <- nrow(p.temp.s)
    p.TS.stat <- p.temp.s %>%
        summarise(
            uid.n = n_distinct(uid),
            uid.ratio = n_distinct(uid) / uid.total
        ) %>%
        arrange(desc(uid.n))

    if (TS.n >= 1) {
        # exact number of TOP species
        TS.total <- p.temp.s %>% n_distinct(sym(species))
        if (TS.n > TS.total) {
            stop(paste0("Can't get more TOP species than total number: ", TS.total))
        } else {
            p.TS.unique <- p.TS.stat %>%
                arrange(desc(uid.n)) %>%
                slice_head(n = TS.n)
        }
    } else {
        p.TS.unique <- p.TS.stat %>% filter(uid.ratio <= TS.n)
    }

    p.TS.unique <- unique(as.vector(unlist(p.TS.unique)))
    # mark top X species
    p %<>% mutate(sym(paste0(prefix, "TS")) := ifelse(sym(species) %in% p.TS.unique, 1, 0))

    # # # # # # # # # #
    # ssos - add species/version columns (0/1)
    # # # # # # # # # #

    ssc <- c("TGOB", "TS", "TO")
    species.unique <- unique(as.vector(unlist(p %>% dplyr::select(sym(species)))))

    # versions
    for (sscn in ssc) {
        ssos.temp.v <- as_tibble(p) %>% dplyr::select(-geometry)
        if (sscn != "TGOB") {
            ssos.temp.v %<>% filter(sym(sscn) == 1)
        }

        # species
        for (sp in species.unique) {
            ssos.temp.v.sp <- ssos.temp.v %>% filter(sym(species) == sp)

            if (nrow(ssos.temp.v.sp) < 1) {
                # zero species occurrences in ssos
                next
            }

            # select observers with at least one occ of selected species
            observers.unique <- unique(as.vector(unlist(ssos.temp.v.sp %>% dplyr::select(sym(observers)))))
            ssos.temp <- ssos.temp.v %>% filter(sym(observers) %in% observers.unique)

            # count occs per species and observers
            ssos.temp.stat <- ssos.temp %>%
                ungroup() %>%
                group_by(sym(species), sym(observers)) %>%
                summarise(uid.n = n_distinct(uid))

            # sum occs per observers
            ssos.temp.stat.observers <- ssos.temp.stat %>%
                ungroup() %>%
                group_by(sym(observers)) %>%
                summarise(observers.n = sum(uid.n))

            # sum occs of selected species per observers
            ssos.temp.stat.species <- ssos.temp.stat %>%
                ungroup() %>%
                filter(sym(species) == sp) %>%
                group_by(sym(observers)) %>%
                summarise(species.n = sum(uid.n))

            # ratio (selected species occs per total occs) per observers
            ssos.temp.ratio <- sssos.temp.stat.observers %>%
                left_join(ssos.temp.stat.species, by = sym(observers)) %>%
                mutate(ratio = species.n / observers.n)

            # remove "outlier" observers with suspicious ratio (lower than centile)
            ssos.temp.ratio.quantile <- unname(quantile(ssos.temp.ratio$ratio, probs = c(0.01)))
            observers.unique <- unique(as.vector(unlist(ssos.temp.ratio %>% filter(ratio > ssos.temp.ratio.quantile[1]) %>% dplyr::select(sym(observers)))))
            p %<>% mutate(sym(paste0(prefix, "ssos", sscn, "_", sp)) := ifelse(sym(observers) %in% observers.unique, 1, 0))

            # IQR version
            limitIQR <- unname(quantile(ssos.temp.ratio$ratio, 0.25)) - (1.5 * IQR(ssos.temp.ratio$ratio))
            observers.unique <- unique(as.vector(unlist(ssos.temp.ratio %>% filter(ratio > limitIQR) %>% dplyr::select(sym(observers)))))
            p %<>% mutate(sym(paste0(prefix, "ssosIqr", sscn, "_", sp)) := ifelse(sym(observers) %in% observers.unique, 1, 0))
        }
    }

    return(p)
}