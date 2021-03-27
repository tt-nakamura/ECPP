NTL = -lntl -lgmp -lm -L/usr/local/lib
ECPP = ECPP.o EC_p.o ZZlib.o rho.o
AM = AtkinMorain.o cornacchia.o IsSquare.o ClassPoly.o CC.o
GK = GoldKil.o EC_pE.o schoof.o
APR = apr.o Cycl.o

test: test.o $(ECPP) $(AM) $(GK)
	g++ test.o $(ECPP) $(GK) $(AM) $(NTL)
fig1: fig1.o $(ECPP) $(AM)
	g++ fig1.o $(ECPP) $(AM) $(NTL)
fig2: fig2.o $(ECPP) $(AM)
	g++ fig2.o $(ECPP) $(AM) $(NTL)
fig3: fig3.o $(ECPP) $(AM) $(GK) $(APR)
	g++ fig3.o $(ECPP) $(AM) $(GK) $(APR) $(NTL)
table1: table1.o ClassNum.o ZZlib.o
	g++ table1.o ClassNum.o ZZlib.o $(NTL)
