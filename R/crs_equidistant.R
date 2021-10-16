#' Create custom equidistant CRS
#'
#' @param sf_obj
#'
#' @return
#' @export
#'
crs_equidistant <- function(sf_obj) {
  wkt_string <-
    'PROJCS["World_Azimuthal_Equidistant",
       GEOGCS["GCS_WGS_1984",
              DATUM["WGS_1984",
              SPHEROID["WGS_1984",6378137,298.257223563]],
              PRIMEM["Greenwich",0],
              UNIT["Degree",0.017453292519943295]],
       PROJECTION["Azimuthal_Equidistant"],
       PARAMETER["False_Easting",0],
       PARAMETER["False_Northing",0],
       PARAMETER["Central_Meridian",{lon_0}],
       PARAMETER["Latitude_Of_Origin",{lat_0}],
       UNIT["Kilometer",1]]'

  central_coords <- sf_obj %>%
    group_by(id) %>% summarize() %>%
    sf::st_centroid() %>%
    sf::st_coordinates() %>% as_tibble()

  lon_0 <- central_coords$X
  lat_0 <- central_coords$Y

  sf::st_crs(glue::glue(wkt_string))
}
