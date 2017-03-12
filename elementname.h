#pragma once
#include<map>

//used as element name storage

#define MAX_ELEMENT_SUPPORT	16
#define MAX_WEIGHT_SUPPORT	64

#define UNKNOWN_WEIGHT	"X"

const string ElementName[MAX_ELEMENT_SUPPORT] = {
	"H","He","Li","Be","B", "c","n","o","F","Ne", "Na","Mg","Al","Si","p", "s"
};
const string ElementNameL[MAX_ELEMENT_SUPPORT] = {
	"H","[He]","[Li]","[Be]","B", "C","N","O","F","[Ne]", "[Na]","[Mg]","[Al]","[Si]","P", "S"
};

const string WeightName[MAX_WEIGHT_SUPPORT] = {
	"", "h","D","T","He","X", "Li","Li","Be","X","X", "B","c","X","n","X", "o","X","X","F","Ne",
	"X","X","Na","Mg","X", "X","Al","Si","X","X", "p","s"
};

const map<int, string> WeightNameL = {
	{1, "H"}, {6, "[Li]"}, {8, "[Be]"}, {12, "C"}, {14, "N"}, {16, "O"}, {19, "F"},
	{23, "[Na]"}, {24, "[Mg]"}, {27, "[Al]"}, {28, "[Si]"}, {31, "P"}, {32, "S"}, {35, "[Cl]"}
};

const map<string, int> NameLWeight = {
	{"H", 1}, {"[Li]", 6}, {"[Be]", 8}, {"B", 11}, {"C", 12}, {"N", 14}, {"O", 16}, {"F", 19},
	{"[Na]", 23}, {"[Mg]", 24}, {"[Al]", 27}, {"[Si]", 28}, {"P", 31}, {"S", 32}, {"[Cl]", 35}
};

const string ShortNames = "HBCNOFPS";

const map<int, int> DefaultHydrogen = {
	{1, 1}, {3, 1}, {8, 2}, {11, 3}, {12, 4}, {14, 3}, {16, 2}, {19, 1},
	{23, 0}, {24, 0}, {27, 3}, {28, 4}, {31, 3}, {32, 2}, {35, 1}
};