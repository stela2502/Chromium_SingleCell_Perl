library(Matrix)
library( RSQLite )

SQLite_ExpressionSummary <- function (fname ) {

	dbh <- dbConnect(SQLite(),dbname=fname )
	sth <- dbSendQuery(dbh, paste(  
			"SELECT gene_id , avg( value), count(value), gname" ,
"from  datavalues left join genes on gene_id = genes.id",
"where sample_id IN (select id from samples where sname not like '%spliced%')",  
"GROUP by gene_id"
		)
	)
	ret <- dbFetch(sth)
	ret
}

SQLite_SampleSummary <- function (fname ) {
	dbh <- dbConnect(SQLite(),dbname=fname )
	sth <- dbSendQuery(dbh, paste(  
					"SELECT sample_id , sum(value) as reads, count(value) as count, sname" ,
					"from  datavalues left join samples on sample_id = samples.id",
					#"where sample_id IN (select id from samples where sname not like '%spliced%')",  
					"GROUP by sample_id"
			)
	)
	ret <- dbFetch(sth)
	ret
}


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

SQLite_2_matrix <- function ( fname, useS=NULL, useG=NULL ) {
	dbh <- dbConnect(SQLite(),dbname=fname )
	ret <- matrix( 
			0, 
			nrow=as.numeric( fetch_first( dbSendQuery(dbh,  "select max(id) from genes" ) ) ),
			ncol=as.numeric( fetch_first( dbSendQuery(dbh,  "select max(id) from samples" ) ) )
	)
	q = "select gene_id, sample_id, value from datavalues where sample_id = :x"

	if ( ! is.null(useG) ) {
		q = paste( q, "and gene_id  IN ( ", paste( collapse=", ", useG),")" )
	}
	sth <- dbSendQuery(dbh, q )
	if(is.null(useS)){
		useS = 1:ncol(ret)
	}
	for ( i in  useS ) {
		dbBind(sth, param = list(x = i))
		t <- dbFetch(sth)
		ret[t$gene_id, i] <- t$value
		if ( i %% 100 == 0 ) {
			print ( paste( "done with sample ",i, "(",nrow(t)," gene entries )"))
		}
	}
	dbClearResult(sth)
	rm(sth)
	q <- "select gname from genes"
	sth <- dbSendQuery(dbh,q)
	rownames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	rm(sth)
	
	sth <- dbSendQuery(dbh,"select sname from samples")
	colnames(ret) <- as.character(t(dbFetch(sth)))
	dbClearResult(sth)
	dbDisconnect(dbh)
	## now tailor the dataset to the with data columns and rows
	ret <- ret[ -which(apply( ret,1,sum) == 0), -which(apply( ret,2,sum) == 0)]
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