SUBDIRS = minhash containmenthash

subdirs:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir;  \
	done
clean: ; rm -rfv ./minhash/bin/minhash* ./minhash/src/*.o ./minhash/src/*.gch ./containmenthash/bin/containmenthash* ./containmenthash/src/*.o ./containmenthash/src/*.gch
