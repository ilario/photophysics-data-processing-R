cedir="ce"
tpvdir="tpv"
tpcdir="tpc"

fromCeToTable(cedir=cedir)
ce(cedir=cedir)
ceFromOutputToGraph(cedir=cedir)
fromTpvTpcToTable(dir=tpvdir)
tpv(tpvdir=tpvdir)
tpvFromOutputToGraph(tpvdir=tpvdir)
tpvCeFromOutputToGraph(cedir=cedir, tpvdir=tpvdir, printcefit=FALSE)
fromTpvTpcToTable(dir=tpcdir)
tpc(tpcdir=tpcdir)
dcFromOutputToGraph(tpvdir=tpvdir, tpcdir=tpcdir)
ceDc(cedir=cedir, tpvdir=tpvdir, tpcdir=tpcdir)
tpcVsTpvVsCe(cedir=cedir, tpvdir=tpvdir, tpcdir=tpcdir)
