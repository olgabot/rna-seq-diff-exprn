#!/usr/bin/Rscript

make.pinkogram = function(){
	# Makes a high-resolution blue-red (blue=low, red=high, can think of it 
	# as cold and hot) color map to represent values in heatmaps. 
	# An alternative to red-green heatmap colors, which are unreadable for
	# those that are color-blind.
	# Author: Olga Botvinnik
	mycol = vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, 
		maxColorValue=255)
	mycol = rev(mycol)
	return(mycol)
}