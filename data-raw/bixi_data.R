library(data.table)

bixi_spatial_features <- fread('data-raw/bixi_spatial_features.csv', header = TRUE, encoding = 'UTF-8')
bixi_spatial_locations <- fread('data-raw/bixi_spatial_locations.csv', header = TRUE, encoding = 'UTF-8')
bixi_station_departures <- fread('data-raw/bixi_station_departures.csv', header = TRUE, encoding = 'UTF-8')
bixi_temporal_features <- fread('data-raw/bixi_temporal_features.csv', header = TRUE, encoding = 'UTF-8')
bixi_temporal_locations <- fread('data-raw/bixi_temporal_locations.csv', header = TRUE, encoding = 'UTF-8')

setkey(bixi_spatial_features, location)
setkey(bixi_spatial_locations, location)
setkey(bixi_station_departures, location)
setkey(bixi_temporal_features, time)
setkey(bixi_temporal_locations, time)

usethis::use_data(
    bixi_spatial_features,
    bixi_spatial_locations,
    bixi_station_departures,
    bixi_temporal_features,
    bixi_temporal_locations,
    overwrite = TRUE
)
