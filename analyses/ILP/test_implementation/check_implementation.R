library(raster)


# Gurobi version with R ---------------------------------------------------

test_R_gurobi651 <- raster("analyses/ILP/test_implementation/test_R_gurobi651.tif")
test_R_gurobi651_include <- raster("analyses/ILP/test_implementation/test_R_gurobi651_include.tif")
test_R_gurobi700 <- raster("analyses/ILP/test_implementation/test_R_gurobi700.tif")

values_R_gurobi651 <- getValues(test_R_gurobi651)
values_R_gurobi651[is.na(values_R_gurobi651)] <- -1

values_R_gurobi651_include <- getValues(test_R_gurobi651_include)
values_R_gurobi651_include[is.na(values_R_gurobi651_include)] <- -1

values_R_gurobi700 <- getValues(test_R_gurobi700)
values_R_gurobi700[is.na(values_R_gurobi700)] <- -1

any(!(values_R_gurobi651 == values_R_gurobi700))
any(!(values_R_gurobi651 == values_R_gurobi651_include))

R_top10 <- test_R_gurobi651 >= 0.9
table(getValues(R_top10))
sum(!is.na(getValues(R_top10)))

# Python implementation ---------------------------------------------------

test_python_gurobi651 <- raster("analyses/ILP/test_implementation/test_python.tif")

python_top10 <- test_python_gurobi651 >= 0.9
table(getValues(python_top10))
sum(!is.na(getValues(python_top10)))
