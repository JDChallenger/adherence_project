para:
	g++ WithinHostModel.cpp
	./a.out
	python plot_parasitaemia.py

para_with_PKPD:
	g++ WithinHostModel_with_PKPD.cpp RK4.cpp
	./a.out
	python plot_parasitaemia_with_PKPD.py

tar:
	cd ..; tar -zcvf within_host_model.tar.gz WithinHostModel