library(Matrix)
library( RSQLite )

SQLite_2_sparseMatrix <- function ( fname ) {
	dbh <- dbConnect(SQLite(),dbname=fname )
	ret <- Matrix( 
			0, 
			nrow=fetch_first( dbSendQuery(dbh,  "select max(id) from genes" ) ),
			ncol=fetch_first( dbSendQuery(dbh,  "select max(id) from samples" ) ),
			sparse = TRUE
	)
	
	sth <- dbSendQuery(dbh,  "select gene_id, sample_id, value from datavalues where sample_id = :x" )
	
	for ( i in 1:ncol(ret)  ) {
		dbBind(sth, param = list(x = i))
		t <- dbFetch(sth)
		if ( nrow(t) > 0 ) {
			ret[t$gene_id, i] <- t$value
		}
		print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
	}
	dbClearResult(sth)
	rm(sth)
	sth <- dbSendQuery(dbh,"select gname from genes")
	rownames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	rm(sth)
	
	sth <- dbSendQuery(dbh,"select sname from samples")
	colnames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	dbDisconnect(dbh)
	ret
}

SQLite_2_matrix <- function ( fname ) {
	dbh <- dbConnect(SQLite(),dbname=fname )
	ret <- matrix( 
			0, 
			nrow=fetch_first( dbSendQuery(dbh,  "select max(id) from genes" ) ),
			ncol=fetch_first( dbSendQuery(dbh,  "select max(id) from samples" ) )
	)
	
	sth <- dbSendQuery(dbh,  "select gene_id, sample_id, value from datavalues where sample_id = :x" )
	
	for ( i in 1:ncol(ret)  ) {
		dbBind(sth, param = list(x = i))
		t <- dbFetch(sth)
		
		ret[t$gene_id, i] <- t$value
		print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
	}
	dbClearResult(sth)
	rm(sth)
	sth <- dbSendQuery(dbh,"select gname from genes")
	rownames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	rm(sth)
	
	sth <- dbSendQuery(dbh,"select sname from samples")
	colnames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	dbDisconnect(dbh)
	ret
}

fetch_first <- function ( sth ) {
	ret <- fetch(sth)
	if (dbHasCompleted(sth)) {
		dbClearResult(sth)
		rm(sth)
	}
	else {
		warn ("query was not finished")
		dbClearResult(sth)
		rm(sth)
	}
	ret
}