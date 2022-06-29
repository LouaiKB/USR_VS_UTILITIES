CC=gcc

all: bin/embed

bin/embed: obj/REAL_preparation.o obj/io_service_pool.o obj/safe_counter.o 
	${CC} -o $@ $^ -L ${RDKIT_ROOT}/lib -lRDKitMolAlign -lRDKitDistGeomHelpers -lRDKitForceFieldHelpers -lRDKitSubstructMatch -lRDKitFileParsers -lRDKitSmilesParse -lRDKitDepictor -lRDKitGraphMol -lRDKitDistGeometry -lRDKitRDGeneral -lRDKitDescriptors -L ${BOOST_ROOT}/lib -lboost_date_time -lboost_thread -lboost_system -lboost_filesystem -lpthread -lm -lstdc++ 

obj/%.o: src/%.cpp
	${CC} -o $@ $^ -c -std=c++17 -O2 -Wall -I ${BOOST_ROOT}/include -I ${BOOST_ROOT}/include

clean:
	rm -f bin/* obj/*
