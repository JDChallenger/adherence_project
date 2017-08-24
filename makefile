para:
	g++ WithinHostModel.cpp
	./a.out
	python plot_parasitaemia.py


tar:
	cd ..; tar -zcvf within_host_model.tar.gz WithinHostModel