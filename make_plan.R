library("drake")
requireNamespace("visNetwork")

r_make(source = "R/drake_plan.R")

r_vis_drake_graph(source = "R/drake_plan.R")

browseURL("kustkarnor.html")
