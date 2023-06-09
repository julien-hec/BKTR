library(data.table)

bixi_spatial_features <- fread('data-raw/bixi_spatial_features.csv', header = TRUE)
bixi_spatial_locations <- fread('data-raw/bixi_spatial_locations.csv', header = TRUE)
bixi_station_departures <- fread('data-raw/bixi_station_departures.csv', header = TRUE)
bixi_temporal_features <- fread('data-raw/bixi_temporal_features.csv', header = TRUE)
bixi_temporal_locations <- fread('data-raw/bixi_temporal_locations.csv', header = TRUE)

setindex(bixi_spatial_features, location)
setindex(bixi_spatial_locations, location)
setindex(bixi_station_departures, location)
setindex(bixi_temporal_features, time)
setindex(bixi_temporal_locations, time)

usethis::use_data(
    bixi_spatial_features,
    bixi_spatial_locations,
    bixi_station_departures,
    bixi_temporal_features,
    bixi_temporal_locations,
    overwrite = TRUE
)
