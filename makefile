bin/bicg : bicg.cpp
	g++ bicg.cpp -o bin/bicg

bin/Front_Bicg : Front_Bicg.c
	gcc Front_Bicg.c -o bin/Front_Bicg