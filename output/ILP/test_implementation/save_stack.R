library(raster)

for (i in 1:length(mcp_ilp)) {
  path <- "analyses/ILP/test_implementation/test_R"
  writeRaster(mcp_ilp[[i]], file.path(path,
                                      paste0(names(mcp_ilp[[i]]), ".tif")))
}
