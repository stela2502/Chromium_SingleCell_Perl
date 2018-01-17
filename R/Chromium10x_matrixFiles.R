to_sparseMatrix <- function( path ){
	if ( ! dir.exists( path ) ) {
		stop( "please - I need an existsing path at start up!")
	}
	ok = 1
	for ( f in c( 'matrix.mtx', 'genes.tsv', 'barcodes.tsv')) {
		if (! file.exists( file.path( path, f) )  ){
			ok = 0
			print ( paste("The required file", f, "is missing in the path") )
		}
	}
	mat = readMM (file.path( path, 'matrix.mtx') )
	genes <- read.delim(file.path( path, 'genes.tsv') )
	cells <-  scan(file.path( path,'barcodes.tsv'), what= character())
	colnames(mat) <- cells
	rownames(mat) = as.vector(genes[,1])
	
}

readMM <- function( filepath ) {
	if (!requireNamespace("stringr", quietly = TRUE)) {
		stop("stringr needed for this function to work. Please install it.",
				call. = FALSE)
	}
	if ( ! dir.exists( filepath ) )  {
		stop( paste( "file", filepath, "not readable") )
	}
	con = file(filepath, "r")
	rows = cols = entries = NULL
	while ( TRUE ) {
		line = readLines(con, n = 1)
		if ( length(line) == 0 ) {
			break
		}
		if ( substr( line, 1,1) == '%'  ) {
			next ## skip comments
		}
		line = unlist(str_split( str, '\\s+'))
		if ( is.null(rows) ) {
			rows = line[1]
			ret <- Matrix( 
					0, 
					nrow=line[1],
					ncol=line[2],
					sparse = TRUE
			)
		}
		ret[line[1], line[2]] = line[3]
	}
	
	close(con)
	ret
}