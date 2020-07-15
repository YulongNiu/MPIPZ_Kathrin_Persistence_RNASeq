datacache <- new.env(hash=TRUE, parent=emptyenv())

org.Ljaponicus.eg <- function() showQCData("org.Ljaponicus.eg", datacache)
org.Ljaponicus.eg_dbconn <- function() dbconn(datacache)
org.Ljaponicus.eg_dbfile <- function() dbfile(datacache)
org.Ljaponicus.eg_dbschema <- function(file="", show.indices=FALSE) dbschema(datacache, file=file, show.indices=show.indices)
org.Ljaponicus.eg_dbInfo <- function() dbInfo(datacache)

org.Ljaponicus.egORGANISM <- "Lotus japonicus"

.onLoad <- function(libname, pkgname)
{
    ## Connect to the SQLite DB
    dbfile <- system.file("extdata", "org.Ljaponicus.eg.sqlite", package=pkgname, lib.loc=libname)
    assign("dbfile", dbfile, envir=datacache)
    dbconn <- dbFileConnect(dbfile)
    assign("dbconn", dbconn, envir=datacache)

    ## Create the OrgDb object
    sPkgname <- sub(".db$","",pkgname)
    db <- loadDb(system.file("extdata", paste(sPkgname,
      ".sqlite",sep=""), package=pkgname, lib.loc=libname),
                   packageName=pkgname)    
    dbNewname <- AnnotationDbi:::dbObjectName(pkgname,"OrgDb")
    ns <- asNamespace(pkgname)
    assign(dbNewname, db, envir=ns)
    namespaceExport(ns, dbNewname)
        
    packageStartupMessage(AnnotationDbi:::annoStartupMessages("org.Ljaponicus.eg.db"))
}

.onUnload <- function(libpath)
{
    dbFileDisconnect(org.Ljaponicus.eg_dbconn())
}

