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
print("Select a DIRECTORY CONTAINING ALL PHOTOPHYSICS DATA OF ONE DEVICE")
datadir=tk_choose.dir(suggestedDataDir,"Select a DIRECTORY CONTAINING ALL PHOTOPHYSICS DATA OF ONE DEVICE")
setwd(datadir)

cedir = file.path(datadir, "ce")
if(!dir.exists(cedir)){
  print("Select CHARGE EXTRACTION Directory")
  cedir=tk_choose.dir(datadir,"Select CHARGE EXTRACTION Directory")
}

tpcdir = file.path(datadir, "tpc")
if(!dir.exists(tpcdir)){
  print("Select TRANSIENT PHOTO CURRENT Directory")
  tpcdir=tk_choose.dir(dirname(cedir),"Select TRANSIENT PHOTO CURRENT Directory")
}

tpvdir = file.path(datadir, "tpv")
if(!dir.exists(tpvdir)){
  print("Select TRANSIENT PHOTO VOLTAGE Directory")
  tpvdir=tk_choose.dir(dirname(cedir),"Select TRANSIENT PHOTO VOLTAGE Directory")
}

source(file.path(scriptsdir,"limits_for_graphics.R"))

source(file.path(scriptsdir,"from_ce_to_table.R"))
source(file.path(scriptsdir,"ce-integrateExp.R"))
source(file.path(scriptsdir,"ce-from_output_to_graph.R"))
source(file.path(scriptsdir,"from_tpv_tpc_to_table.R"))
source(file.path(scriptsdir,"tpv.R"))
source(file.path(scriptsdir,"tpv-from_output_to_graph.R"))
source(file.path(scriptsdir,"tpvce-from_output_to_graph.R"))
source(file.path(scriptsdir,"tpc.R"))
source(file.path(scriptsdir,"dc-from_output_to_graph.R"))
source(file.path(scriptsdir,"ce-with_limits.R"))
source(file.path(scriptsdir,"cedc.R"))
source(file.path(scriptsdir,"tpc-vs-tpv-vs-ce.R"))
source(file.path(scriptsdir,"jrec-ce.R"))
source(file.path(scriptsdir,"jrec-dc.R"))
source(file.path(scriptsdir,"tpv-from_output_to_graph-with_limits.R"))

fromCeToTable(cedir=cedir)
ceIntegrateExp(cedir=cedir)
ceFromOutputToGraph(cedir=cedir)
fromTpvTpcToTable(dir=tpvdir)
tpv(tpvdir=tpvdir)
tpvFromOutputToGraph(tpvdir=tpvdir)
tpvCeFromOutputToGraph(cedir=cedir, tpvdir=tpvdir, printcefit=FALSE)
fromTpvTpcToTable(dir=tpcdir)
tpc(tpcdir=tpcdir)
dcFromOutputToGraph(tpvdir=tpvdir, tpcdir=tpcdir)
ceWithLimits(cedir=cedir, dcdir=getwd())
ceDc(cedir=cedir, tpvdir=tpvdir, tpcdir=tpcdir)
tpcVsTpvVsCe(cedir=cedir, tpvdir=tpvdir, tpcdir=tpcdir)
jrecCe(cedir=cedir, tpvdir=tpvdir)
jrecDc(tpvdir=tpvdir)
tpvFromOutputToGraphWithLimits(tpvdir=tpvdir, dcdir=getwd())
