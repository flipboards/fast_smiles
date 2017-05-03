#include"smiles.h"

int main() {
	printf("test for smiles...\n");
	smiles m;
	m.read_smiles("[C1H2][CH2][CH2]1");
	/*m.push_atom(16, 1);
	m.push_atom(16, 1);
	m.push_atom(16, 1);
	m.connect_atom(0, 1, 2);
	m.connect_atom(1, 2, 2);*/

	//m.push_atom(1);
	m.canonicalize();
	m.compile(true, true);
	printf("%s", m.to_string().c_str());
	return 0;
}