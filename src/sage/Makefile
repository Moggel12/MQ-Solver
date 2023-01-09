all: run

fes_rec.py: fes_rec.sage
	sage --preparse fes_rec.sage
	mv fes_rec.sage.py fes_rec.py

fes.py: fes.sage
	sage --preparse fes.sage
	mv fes.sage.py fes.py

mob_new.py: mob_new.sage
	sage --preparse mob_new.sage
	mv mob_new.sage.py mob_new.py

utils.py: utils.sage
	sage --preparse utils.sage
	mv utils.sage.py utils.py

run: fes_rec.py fes.py mob_new.py utils.py
	sage dinur.sage
	
.PHONY: clean

clean:
	rm -f fes.py fes_rec.py mob_new.py utils.py dinur.sage.py
