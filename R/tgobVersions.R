#' @title Subsample versions (TO, TS, TSAO, TOsso, TSsso, TSAOsso, TGOBsso) from TGOB
#'
#' @param p sf/sfc (POINT/MULTIPOINT): point occurrences (TGOB)
#' @param r RasterLayer: Template raster. If provided, all stats are calculated per raster pixel (cell/square).
#' @param species character: column name with species name
#' @param observers character: column name with observers' names (expected: "John Doe" / "Doe J." / "John Doe, Jane Maria Moe" / "Doe J., Moe J. M." / ...)
#' @param TS.n integer or numeric: int or ratio (< 1, \code{cumsum} \emph{p}) to select for TOP species
#' @param TO.n integer or numeric: int or ratio (< 1, \code{cumsum} \emph{p}) to select for TOP observers
#' @param TSAO.n integer or numeric: int or quantile (< 1, upper quartile 0.75, ...) to select for TOP species by Area Overlap (overlaping occupied cells with focus species). Can be done only with inserted \emph{r}!
#' @param TSAO.min numeric (0-1): minimum ratio overlap to threshold before \emph{TSAO.n} is applied
#' @param prefix character: prefix for newly added columns to \emph{p} with binary (0/1) sign that links affiliation to version (to select it easily later)
#' @param observersRemoveSingleName logical: remove single word (name) observers
#' @param observersRemoveSingleOccurrence integer: remove observers with single (or specified number) species occurrence from sso. Can be too restrictive for very rare species, where observers have small chance to re-observe. On the other hand, can be useful to restrict very active observers that record some species only purposefully.
#' @param quantileObserversThreshold numeric: (0-1) (centile 0.01, decile 0.10, ...) threshold to remove individual lower outlier observers with "suspicious" (purposefull) pattern, calculated as \emph{observational ratio} (Σ\emph{focus species} / Σ\emph{rest of species}) compared to overall observers \emph{observational ratio}. In other words: comparing (per species, \emph{sso}?) each observers' profile (sums of species occurerrences) to overall observers' profile and then remove outliers.
#' @param badWordsSpecies character vector: if species name containing unwanted strings, remove such rows
#' @param badWordsObservers character vector: if observers name containing unwanted strings, remove such rows
#' @param crs integer: Force crs.
#'
#' @return named list: \emph{p}: sf (POINT/MULTIPOINT): inserted \emph{p} with new (0/1) columns: TO, TS, TSAO ssoX_\emph{species}; \emph{report}: TO+TS+TSAO stats used to threshold selected top O+S
#'
#' @export

tgobVersions <- function(p, r = NA, species = "species", observers = NA, TS.n = 0.2, TO.n = 0.2, TSAO.n = 0.8, TSAO.min = 0.5, prefix = "nc_", observersRemoveSingleName = TRUE, observersRemoveSingleOccurrence = 0, quantileObserversThreshold = 0, badWordsSpecies = NA, badWordsObservers = NA, crs = NA) {

    # not change original data, just add new prefixed columns with weights or binary sign

    minMaxNormalize <- function(v) {
        v.min <- min(v)
        v.max <- max(v)
        return((v - v.min) / (v.max - v.min))
    }

    badWords <- "_badWords"
    uid <- paste0(prefix, "uid")

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

    p %<>% ungroup() %>% mutate(!!uid := row_number())

    # # # # # # # # # #
    # prepare sso  + prefilter badWords
    # # # # # # # # # #

    if (!is.na(observers)) {
        # re-save observers' names to prefixed column and re-assign original variable name
        observers.temp <- paste0(prefix, observers)
        p %<>% mutate(!!observers.temp := !!sym(observers))
        observers <- observers.temp

        #
        # vše níže v nastavení - nechat základní vhodné defaultní a další volitelné,,, udělat nějak hromadně vektorem, co chci odstranit
        #
        # replace accents
        p %<>% mutate(!!observers := stringi::stri_trans_general(str = !!sym(observers), id = "Any-ASCII"))

        # lowercase
        p %<>% mutate(!!observers := stringi::stri_trans_general(str = !!sym(observers), id = "Any-Lower"))

        # make two words from selected character (not only dot TODO) separated strings
        p %<>% mutate(!!observers := str_replace(!!sym(observers), "\\.", "\\. "))

        # trim spaces and internal multi spaces
        p %<>% mutate(!!observers := stringr::str_squish(!!sym(observers)))

        # replace non alphanumeric+space chars - možnost zadat přímo string do reguláru!!!
        p %<>% mutate(!!observers := str_replace_all(!!sym(observers), "[^[:alnum:]| ]", "_"))


        # mark single word observers (nesmysl, pokud jsou tam unikátní stringy, třeba přezdívky)
        p %<>% mutate(!!paste0(observers, "_single") := ifelse(stringi::stri_count_fixed(!!sym(observers), " ") == 0, 1, 0))
        if (observersRemoveSingleName) {
            # remove single word observers
            p %<>% filter(!!sym(paste0(observers, "_single")) == 0)
        }

        if (!is.na(badWordsObservers)) {
            # mark bad words observers
            p %<>% mutate(!!paste0(observers, badWords) := ifelse(str_detect(!!sym(observers), paste(badWordsObservers, collapse = "|")), 1, 0))
            p %<>% filter(!!sym(paste0(observers, badWords)) == 0)
        } else {
            p %<>% mutate(!!paste0(observers, badWords) := 0)
        }
    }

    if (!is.na(badWordsSpecies)) {
        # mark bad words species
        p %<>% mutate(!!paste0(prefix, species, badWords) := ifelse(str_detect(!!sym(species), paste(badWordsSpecies, collapse = "|")), 1, 0))
        p %<>% filter(!!sym(paste0(prefix, species, badWords)) == 0)
    } else {
        p %<>% mutate(!!paste0(species, badWords) := 0)
    }



    #  dát možnost volby nechat výstup dofiltrovaný, nebo si ho manuálně ručně dofiltrovat sám až potom
    # - ale nezahrnovat je do výpočtu statistik!!!

    # nedělal jsem pro celkové TGOB zde normalizaci per pixel, přestože u TS a TO jsem ji udělal!! - jo, ale statistiky jsem zatím dělal jen na d TO a TS, je to OK...
    # když ale zadám raster, tak ty výpočty chci mít normalizované!!! Ale původní dataset bych neměl měnit co do počtu řádku, jen nevhodné bu%nky dostanou 0
    # prs pixel cellNumber ale můžu udělat vždy až dodatečně níže v kódu, není nutné mít tady

    # místo rasteru můžu použít i buffer čistě pro statistíky nad presencemi - pak budu ale budu muset řešit prostorový překryv přímo prostorovou funkcí - ne jen pomocí cellNumber
    if (is(r, "RasterLayer")) {
        #
        # per species, (observer) and pixel
        #
        cellNumber <- paste0(prefix, "cellNumber")
        v.temp <- raster::extract(r, sf::st_coordinates(p), cellnumbers = TRUE)
        p %<>% mutate(!!cellNumber := v.temp[, 1])

        occs.na <- nrow(p %>% filter(is.na(!!sym(cellNumber))))

        if (occs.na > 0) {
            message(paste0(occs.na, " occurences removed, no cellnumbers (raster::extract)"))
            p %<>% filter(!is.na(!!sym(cellNumber)))
        }

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(!!sym(observers), !!sym(species), !!sym(cellNumber), !!sym(uid), -geometry)) %>%
                ungroup() %>%
                group_by(!!sym(observers), !!sym(cellNumber), !!sym(species)) %>% # původně bez species, pak to bylo ale jen z hlediska prostorové aktivity, ale bez počtu sesbíraných druhů
                slice_head(n = 1)
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(!!sym(species), !!sym(cellNumber), !!sym(observers), !!sym(uid), -geometry)) %>%
            ungroup() %>%
            group_by(!!sym(species), !!sym(cellNumber), !!sym(observers)) %>% # původně bez observers, pak to bylo ale jen z hlediska prostorového rozšíření, ale bez počtu sesbíraných druhů
            slice_head(n = 1)
    } else {
        #
        # per species and (observer)
        #

        # TO
        if (!is.na(observers)) {
            # per observer stats
            p.temp.o <- as_tibble(p %>% dplyr::select(!!sym(observers), !!sym(uid), -geometry)) %>%
                group_by(!!sym(observers))
        }

        # TS
        # per species stats
        p.temp.s <- as_tibble(p %>% dplyr::select(!!sym(species), !!sym(uid), -geometry)) %>%
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
                uid.n = n_distinct(!!sym(uid)),
                uid.ratio = n_distinct(!!sym(uid)) / uid.total
            ) %>%
            arrange(desc(uid.n))

        if (TO.n >= 1) {
            # exact number of TOP observers
            TO.total <- p.temp.o %>% n_distinct(!!sym(observers))
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

        ### # wTO
        p.TO.stat.w <- p.TO.stat %>%
            mutate(uid.ratio = uid.ratio^2) %>%
            dplyr::select(!!sym(observers), uid.ratio) %>%
            rename(!!paste0(prefix, "TO.w") := "uid.ratio")

        p.TO.stat.w[[paste0(prefix, "TO.w")]] <- minMaxNormalize(p.TO.stat.w[[paste0(prefix, "TO.w")]])

        p %<>% left_join(p.TO.stat.w, by = observers)
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
            uid.n = n_distinct(!!sym(uid)),
            uid.ratio = n_distinct(!!sym(uid)) / uid.total
        ) %>%
        arrange(desc(uid.n))

    if (TS.n >= 1) {
        # exact number of TOP species
        TS.total <- p.temp.s %>% n_distinct(!!sym(species))
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

    ### # wTS
    p.TS.stat.w <- p.TS.stat %>%
        mutate(uid.ratio = uid.ratio^2) %>%
        dplyr::select(!!sym(species), uid.ratio) %>%
        rename(!!paste0(prefix, "TS.w") := "uid.ratio")

    p.TS.stat.w[[paste0(prefix, "TS.w")]] <- minMaxNormalize(p.TS.stat.w[[paste0(prefix, "TS.w")]])

    p %<>% left_join(p.TS.stat.w, by = species)

    # # # # # # # # # #
    # sso - add species/version columns (0/1)
    # # # # # # # # # #


    species.unique <- unique(as.vector(unlist(as_tibble(p) %>% dplyr::select(!!sym(species)))))

    # sso versions
    p.t <- as_tibble(p) %>% dplyr::select(-geometry)
    TSAO.min.species <- c()
    TOAO.min.species <- c()
    TSAO.top <- list()
    TOAO.top <- list()

    # species
    for (sp in species.unique) {
        sso.temp <- NA
        message(sp)
        p.t.sp <- p.t %>% filter(!!sym(species) == sp)

        if (nrow(p.t.sp) < 1) {
            # zero species occurrences in sso
            next
        }

        # select observers with at least one occ of selected species
        observers.unique <- unique(as.vector(unlist(p.t.sp %>% dplyr::select(!!sym(observers)))))
        sso.temp <- p.t %>% filter(!!sym(observers) %in% observers.unique)

        # count occs per species and observers
        sso.temp.stat <- sso.temp %>%
            ungroup() %>%
            group_by(!!sym(species), !!sym(observers)) %>%
            summarise(uid.n = n_distinct(!!sym(uid)))

        # sum occs per observers
        sso.temp.stat.observers <- sso.temp.stat %>%
            ungroup() %>%
            group_by(!!sym(observers)) %>%
            summarise(observers.n = sum(uid.n))

        # sum occs of selected species per observers
        sso.temp.stat.species <- sso.temp.stat %>%
            ungroup() %>%
            filter(!!sym(species) == sp) %>%
            group_by(!!sym(observers)) %>%
            summarise(species.n = sum(uid.n))

        # ratio (selected species occs per total occs) per observers
        # počet (per pixel) pozorování pro daného pozorovatele: species.n (focus druhu) / observers.n (celkem všech druhů) = ratio
        sso.temp.ratio <- sso.temp.stat.observers %>%
            left_join(sso.temp.stat.species, by = as.character(sym(observers))) %>%
            mutate(ratio = species.n / observers.n)

        observersRemoveSingleOccurrence <- as.numeric(observersRemoveSingleOccurrence)
        if (observersRemoveSingleOccurrence > 0) {
            # remove single
            sso.temp.ratio %<>% filter(species.n > observersRemoveSingleOccurrence)
            observers.unique <- unique(as.vector(unlist(sso.temp.ratio %>% dplyr::select(!!sym(observers)))))
        }

        quantileObserversThreshold <- as.numeric(quantileObserversThreshold)
        if (quantileObserversThreshold > 0 & quantileObserversThreshold < 1) {
            # remove quantile (and/or singles before)
            # remove "outlier" observers with suspicious ratio (lower than input quantile)
            sso.temp.ratio.quantile <- unname(stats::quantile(sso.temp.ratio$ratio, probs = quantileObserversThreshold, na.rm = TRUE))
            observers.unique <- unique(as.vector(unlist(sso.temp.ratio %>% filter(ratio > sso.temp.ratio.quantile[1]) %>% dplyr::select(!!sym(observers)))))
        }

        p %<>% mutate(!!paste0(prefix, "sso", "_", sp) := ifelse(!!sym(observers) %in% observers.unique, 1, 0))


        ### # wsso
        sso.temp.ratio.w <- sso.temp.ratio %>%
            dplyr::select(!!sym(observers), ratio) %>%
            rename(!!paste0(prefix, "TGOB.sso.w", "_", sp) := "ratio")

        sso.temp.ratio.w[[paste0(prefix, "TGOB.sso.w", "_", sp)]] <- minMaxNormalize(sso.temp.ratio.w[[paste0(prefix, "TGOB.sso.w", "_", sp)]])

        p %<>% left_join(sso.temp.ratio.w, by = observers)


        ### # wsso - alts
        sso.temp.ratio[["observers.n.n"]] <- minMaxNormalize(sso.temp.ratio$observers.n)
        sso.temp.ratio[["species.n.n"]] <- minMaxNormalize(sso.temp.ratio$species.n)
        sso.temp.ratio[[paste0(prefix, "TGOB.sso.w.3", "_", sp)]] <- minMaxNormalize(1 / 3 * sso.temp.ratio$ratio + 1 / 3 * sso.temp.ratio$observers.n.n + 1 / 3 * sso.temp.ratio$species.n.n)

        sp.w <- sso.temp.ratio
        sp.w %<>% arrange(ratio)
        sp.sort <- sort(sp.w$ratio)
        sp.mean <- mean(sp.w$ratio)
        sp.median <- median(sp.w$ratio)
        sp.sd <- sd(sp.w$ratio)

        # dnorm1 median
        n <- dnorm(1:length(sp.sort), mean = which.min(abs(sp.sort - sp.median)), sd = sd(1:length(sp.sort)))
        sp.w[[paste0(prefix, "TGOB.sso.w.dnorm1.median", "_", sp)]] <- round(minMaxNormalize(n), digits = 5)
        sp.w[[paste0(prefix, "TGOB.sso.w.3.dnorm1.median", "_", sp)]] <- minMaxNormalize(1 / 3 * sp.w[[paste0(prefix, "TGOB.sso.w.dnorm1.median", "_", sp)]] + 1 / 3 * sp.w$observers.n.n + 1 / 3 * sp.w$species.n.n)


        # dnorm1 mean
        n <- dnorm(1:length(sp.sort), mean = which.min(abs(sp.sort - sp.mean)), sd = sd(1:length(sp.sort)))
        sp.w[[paste0(prefix, "TGOB.sso.w.dnorm1.mean", "_", sp)]] <- round(minMaxNormalize(n), digits = 5)
        sp.w[[paste0(prefix, "TGOB.sso.w.3.dnorm1.mean", "_", sp)]] <- minMaxNormalize(1 / 3 * sp.w[[paste0(prefix, "TGOB.sso.w.dnorm1.mean", "_", sp)]] + 1 / 3 * sp.w$observers.n.n + 1 / 3 * sp.w$species.n.n)


        # dnorm2 median
        n <- dnorm(sp.sort, mean = sp.median, sd = sp.sd)
        sp.w[[paste0(prefix, "TGOB.sso.w.dnorm2.median", "_", sp)]] <- round(minMaxNormalize(n), digits = 5)

        sp.w[[paste0(prefix, "TGOB.sso.w.3.dnorm2.median", "_", sp)]] <- minMaxNormalize(1 / 3 * sp.w[[paste0(prefix, "TGOB.sso.w.dnorm2.median", "_", sp)]] + 1 / 3 * sp.w$observers.n.n + 1 / 3 * sp.w$species.n.n)

        # dnorm2 mean
        n <- dnorm(sp.sort, mean = sp.mean, sd = sp.sd)
        sp.w[[paste0(prefix, "TGOB.sso.w.dnorm2.mean", "_", sp)]] <- round(minMaxNormalize(n), digits = 5)

        sp.w[[paste0(prefix, "TGOB.sso.w.3.dnorm2.mean", "_", sp)]] <- minMaxNormalize(1 / 3 * sp.w[[paste0(prefix, "TGOB.sso.w.dnorm2.mean", "_", sp)]] + 1 / 3 * sp.w$observers.n.n + 1 / 3 * sp.w$species.n.n)

        p %<>% left_join(sp.w[, c(1, 7:15)], by = observers) # TODO get cols by name + autors to join

        # # # # # # # # # #
        # TSAO
        # # # # # # # # # #
        if (is(r, "RasterLayer")) {
            ss.total <- length(unique(p.t.sp[[cellNumber]]))

            p.TSAO.stat <- p.temp.s %>%
                ungroup() %>%
                group_by(!!sym(species)) %>%
                filter(!!sym(cellNumber) %in% unique(p.t.sp[[cellNumber]])) %>%
                summarise(
                    cells.shared = n_distinct(!!sym(cellNumber)),
                    cells.shared.ratio = n_distinct(!!sym(cellNumber)) / ss.total
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
            p %<>% mutate(!!paste0(prefix, "TS.AO", "_", sp) := ifelse(!!sym(species) %in% p.TSAO.unique, 1, 0))

            ### # wTSAO
            p.TSAO.stat.w <- p.TSAO.stat %>%
                mutate(cells.shared.ratio = cells.shared.ratio^2) %>%
                dplyr::select(!!sym(species), cells.shared.ratio) %>%
                rename(!!paste0(prefix, "TS.AO.w", "_", sp) := "cells.shared.ratio")

            p.TSAO.stat.w[[paste0(prefix, "TS.AO.w", "_", sp)]] <- minMaxNormalize(p.TSAO.stat.w[[paste0(prefix, "TS.AO.w", "_", sp)]])

            p %<>% left_join(p.TSAO.stat.w, by = species)
        }

        # # # # # # # # # #
        # TOAO
        # # # # # # # # # #
        if (is(r, "RasterLayer")) {
            ### v názvech proměnných mám species místo observers (významově hledám nálezy pozorovbatelů, ne druhů...)
            ss.total <- length(unique(p.t.sp[[cellNumber]]))

            p.TOAO.stat <- p.temp.o %>%
                ungroup() %>%
                group_by(!!sym(observers)) %>%
                filter(!!sym(cellNumber) %in% unique(p.t.sp[[cellNumber]])) %>%
                summarise(
                    cells.shared = n_distinct(!!sym(cellNumber)),
                    cells.shared.ratio = n_distinct(!!sym(cellNumber)) / ss.total
                ) %>%
                arrange(desc(cells.shared))

            # remove focus observers (100 % overlap...)
            # p.TOAO.stat %<>% filter(!!sym(observers) != sp)

            if (nrow(p.TOAO.stat) > 0) {
                if (TSAO.n >= 1) {
                    # exact number of TOP species
                    p.TOAO.unique <- p.TOAO.stat %>% slice_head(n = TSAO.n)
                } else {
                    if (nrow(p.TOAO.stat %>% filter(cells.shared.ratio > TSAO.min)) >= 4) {
                        # left only >50 % overlaping species with more than 3 species
                        p.TOAO.stat.temp <- p.TOAO.stat %>% filter(cells.shared.ratio > TSAO.min)
                        toao.quantile <- NA
                        toao.quantile <- unname(stats::quantile(p.TOAO.stat.temp$cells.shared, probs = TSAO.n, na.rm = TRUE))

                        if (nrow(p.TOAO.stat.temp %>% filter(cells.shared > toao.quantile)) >= 4) {
                            p.TOAO.unique <- p.TOAO.stat.temp %>% filter(cells.shared > toao.quantile)
                        } else {
                            message(paste0("Less than 3 ovelaping observers with >", TSAO.min, " overlap and quantile ", TSAO.n, ". Left 3 top observers"))
                            TOAO.min.species <- c(TOAO.min.species, sp)
                            p.TOAO.unique <- p.TOAO.stat %>% slice_head(n = 4)
                        }
                    } else {
                        message(paste0("Less than 3 ovelaping observers with >", TSAO.min, " overlap. Left 3 top observers"))
                        TOAO.min.species <- c(TOAO.min.species, sp)
                        p.TOAO.unique <- p.TOAO.stat %>% slice_head(n = 4)
                    }
                    TOAO.top[[sp]] <- p.TOAO.unique
                }
            } else {
                message("No overlaping species...")
                next
            }

            p.TOAO.unique <- unique(as.vector(unlist(p.TOAO.unique %>% dplyr::select(!!sym(observers)))))

            # mark top X  AO species
            p %<>% mutate(!!paste0(prefix, "TO.AO", "_", sp) := ifelse(!!sym(observers) %in% p.TOAO.unique, 1, 0))

            ### # wTOAO
            p.TOAO.stat.w <- p.TOAO.stat %>%
                mutate(cells.shared.ratio = cells.shared.ratio^2) %>%
                dplyr::select(!!sym(observers), cells.shared.ratio) %>%
                rename(!!paste0(prefix, "TO.AO.w", "_", sp) := "cells.shared.ratio")

            p.TOAO.stat.w[[paste0(prefix, "TO.AO.w", "_", sp)]] <- minMaxNormalize(p.TOAO.stat.w[[paste0(prefix, "TO.AO.w", "_", sp)]])

            p %<>% left_join(p.TOAO.stat.w, by = observers)
        }
    }

    return(list(
        "p" = p,
        "report" = list(
            "TS" = p.TS.stat, "TS.n" = TS.n,
            "TO" = p.TO.stat, "TO.n" = TO.n,
            "TSAO" = TSAO.top, "TSAO.n" = TSAO.n, "TOAO" = TOAO.top,
            "TSAO.min.species" = TSAO.min.species, "TOAO.min.species" = TOAO.min.species
        )
    ))
}