
CLA=clang++
CXX=g++
CXXFLAGS=-std=c++11 -Ofast -DFINAL_CHECK -DSPECIAL_HP -fpermissive
DEPS=./optimizer.cpp ./optimizer.h ./Utils/reader.h ./Utils/network.h ./Utils/codon.h ./Utils/energy_parameter.h ./Utils/intl11.h ./Utils/intl21.h ./Utils/intl22.h ./Utils/utility_v.h ./Utils/common.h ./Utils/base.h 
BIN=bin/mRNADesign 

mrnadesign: $(DEPS)
	@echo "Compiling" $@ "from" $< "..."
	mkdir -p ./bin
	if $(CLA) $(CXXFLAGS) ./mRNA_design.cpp -o $(BIN); then \
		echo "Compiled with clang++; finished\n"; \
	else \
		if $(CXX) $(CXXFLAGS) ./mRNA_design.cpp -o $(BIN); then \
			echo "Clang++ not available; compiled with g++; finished\n"; \
		else \
			echo "Compilation failed!\n"; \
		fi \
	fi

.PHONY : clean	

clean:
	rm -f $(BIN)
