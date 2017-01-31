#include"smiles.h"

int main() {
	printf("test for smiles...\n");
	smiles m;
	m.read_smiles("[C1H2][CH2][CH2]1");

	m.canonicalize();
	m.compile(true);
	printf("%s", m.to_string().c_str());
	return 0;
}