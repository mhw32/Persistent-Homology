#!/usr/bin/r

readBin <- function(infile) {
	# read the binary file and reshape the
	# matrix into the right dimensions.
	con <- file(infile, "rb")
	dim <- readBin(con, "integer", 2)
	Mat <- matrix( readBin(con, "numeric", prod(dim)), dim[1], dim[2])
	close(con)
	return Mat
}



