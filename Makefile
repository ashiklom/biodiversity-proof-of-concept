.PHONY: drake getadapt

all: drake

drake:
	Rscript -e "drake::r_make()"

getadapt:
	rsync -avz --progress adapt:~/projects/usgs-biodiversity/mstmip-stats.tif data/
