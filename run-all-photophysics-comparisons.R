library(tcltk)

homedir=path.expand('~')#"C:\\TPV_Controls_ORGANICS"

#if(!exists("scriptsdir")){
print("Select R SCRIPTS Directory")
scriptsdir = tk_choose.dir(homedir,"Select R SCRIPTS Directory\n")#"~/software/photophysics-data-processing-R"#"C:\\Users\\iciq\\Desktop\\photophysics-data-processing-R"
#}

print("Select a DIRECTORY CONTAINING ALL DEVICES' PHOTOPHYSICS DATA DIRECTORIES")
datadir=tk_choose.dir(homedir,"Select a DIRECTORY CONTAINING ALL PHOTOPHYSICS DATA OF ONE DEVICE\n")
setwd(datadir)

# cedir = file.path(datadir, "ce")
# if(!dir.exists(cedir)){
#   print("Select CHARGE EXTRACTION Directory")
#   cedir=tk_choose.dir(homedir,"Select CHARGE EXTRACTION Directory\n")
# }
# 
# tpcdir = file.path(datadir, "tpc")
# if(!dir.exists(tpcdir)){
#   print("Select TRANSIENT PHOTO CURRENT Directory")
#   tpcdir=tk_choose.dir(dirname(cedir),"Select TRANSIENT PHOTO CURRENT Directory\n")
# }
# 
# tpvdir = file.path(datadir, "tpv")
# if(!dir.exists(tpvdir)){
#   print("Select TRANSIENT PHOTO VOLTAGE Directory")
#   tpvdir=tk_choose.dir(dirname(cedir),"Select TRANSIENT PHOTO VOLTAGE Directory\n")
# }

source(file.path(scriptsdir,"ce-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"ce-time-many.R"))
source(file.path(scriptsdir,"tpv-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"tpvce-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"dc-from_many_output_to_graph_capacitance.R"))
source(file.path(scriptsdir,"dc-from_many_output_to_graph_charge.R"))
source(file.path(scriptsdir,"tpvdc-from_many_output_to_graph.R"))


