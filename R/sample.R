#' @title Sample from raster or points
#'
#' @param rp RasterLayer or sf/sfc (POINT/MULTIPOINT): Raster or points.
#' @param n integer or numeric: Ratio (< 1) or int to sample.
#' @param prob logical: For raster with probabilities as pixel values only.
#' @param crs integer: Force crs.
#'
#' @references Inspired by: https://github.com/danlwarren/ENMTools/blob/004a4a1e182127900a5f62bc015770479bcd0415/R/check.bg.R#L138-L144
#' @return sf (POINT/MULTIPOINT)
#' 
#' @export

sample <- function(rp, n = 0.1, prob = FALSE, crs = NA) {
    # počítám s tím, že dostanu už unikátní body per pixel! Pak můžu vždy poměrově vybrat část k "ploše"?
    # raster::sampleRandom have some issues not ideal to use in this case?

    if (is(rp, "RasterLayer")) {
        crs <- if (!is.integer(crs)) raster::crs(rp)
        background.points <- as.data.frame(raster::rasterToPoints(rp))
    } else if (is(rp, "sf") && sf::st_geometry_type(rp, by_geometry = FALSE) == "POINT") {
        if (prob != FALSE) {
            stop("rp: sf/sfc (POINT/MULTIPOINT) with prob=FALSE allowed!")
            return(list())
        }
        crs <- if (!is.integer(crs)) sf::st_crs(rp)
        background.points <- as.data.frame(sf::st_coordinates(rp))
        colnames(background.points) <- c("x", "y")
    } else {
        stop("rp: only RasterLayer or sf/sfc (POINT/MULTIPOINT) allowed!")
        return(list())
    }

    if (n >= 1) {
        # exact number of points to sample
        total <- nrow(background.points)
        if (n > total) {
            print(paste0("Can't sample more than total number! ", total, " used, instead of n=", n))
            proportion <- total
        } else {
            proportion <- n
        }
    } else {
        # ratio
        total <- nrow(background.points)
        proportion <- round(total * n)
    }

    probs <- NULL
    if (prob == TRUE) {
        probs <- background.points[, 3]
    }
    inds <- base::sample(1:total,
        size = proportion,
        prob = probs,
        replace = FALSE
    )
    bg.temp <- background.points[inds, 1:2]
    return(bg.temp %>% sf::st_as_sf(coords = c("x", "y"), crs = crs))
}