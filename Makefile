.PHONY: drake

all: drake

drake:
	Rscript -e "drake::r_make()"
