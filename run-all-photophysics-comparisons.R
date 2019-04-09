library(tcltk)

homedir=path.expand('~')#"C:\\TPV_Controls_ORGANICS"

if(exists("scriptsdir")){
  suggestedScriptsDir = scriptsdir
} else {
  suggestedScriptsDir = homedir
}

print("Select photophysics-data-processing.R Directory")
scriptsdir = tk_choose.dir(suggestedScriptsDir,"Select photophysics-data-processing.R Directory")#"~/software/photophysics-data-processing-R"#"C:\\Users\\iciq\\Desktop\\photophysics-data-processing-R"

if(exists("datadir")){
  suggestedDataDir = datadir
} else {
  suggestedDataDir = homedir
}

print("Select a DIRECTORY CONTAINING ALL DEVICES' PHOTOPHYSICS DATA DIRECTORIES")
datadir=tk_choose.dir(suggestedDataDir,"Select a DIRECTORY CONTAINING ALL DEVICES' PHOTOPHYSICS DATA DIRECTORIES")
setwd(datadir)

source(file.path(scriptsdir,"limits_for_graphics.R"))

source(file.path(scriptsdir,"ce-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"ce-time-many.R"))
source(file.path(scriptsdir,"tpv-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"tpvce-from_many_output_to_graph.R"))
source(file.path(scriptsdir,"dc-from_many_output_to_graph_capacitance.R"))
source(file.path(scriptsdir,"dc-from_many_output_to_graph_charge.R"))
source(file.path(scriptsdir,"tpvdc-from_many_output_to_graph.R"))


