para:
	g++ WithinHostModel.cpp
	./a.out
	python plot_parasitaemia.py

#26th November: code has been updated, so that program runs faster
#Original code can still be accessed (see command in comments below)
para_with_PKPD:
	g++ WithinHostModel_with_PKPD_faster.cpp RK4.cpp
	./a.out
	python plot_parasitaemia_with_PKPD.py

#para_with_PKPD_original:
#	g++ WithinHostModel_with_PKPD.cpp RK4.cpp
#	./a.out
#	python plot_parasitaemia_with_PKPD.py

para_with_adherence:
	g++ WithinHostModel_with_adherence.cpp RK4.cpp
	./a.out

tar:
	cd ..; tar -zcvf within_host_model.tar.gz WithinHostModel