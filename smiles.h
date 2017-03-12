/*
A canonical SMILES implementation.
Last modified in 2017, by Yanze Wu.
*/

#pragma once

#ifndef SMILES_H
#define SMILES_H

#include <vector>
#include <list>
#include <map>

#define maxhnum		0x08
#define maxnhcon	0x08
#define maxatmwht	0x80

using namespace std;

class smiles {
public:
	smiles() {}

	//push an atom back with smiles score. do not use directly.
	void _push_back(int score);

	//push a child to an atom without changing score. do not use directly.
	void _push_child(int parent, int child);

	//push an atom back with weight, hydrogennumber (-1 for default hydrogen)
	void push_atom(int weight, int hydrogen=-1);

	//connect a parent with child
	void connect_atom(int atom1, int atom2, bool change_hydrogen=true);

	//number of atoms
	size_t size()const { return atoms.size(); }

	//return the number of atom with certain atom weight
	int element_num(int weight)const;

	bool operator==(const smiles&)const;


	//turning a smiles into a canonical smiles
	void canonicalize();
	
	//Write the string form in symbol
	void compile(bool hydrogen=true);

	//output smiles without check if compiled before.
	string to_string()const { return symbol; }

	//output smiles.
	string to_smiles(bool hydrogen=true);

	//read smiles string. This function may contain error.
	void read_smiles(const string& sms, bool read_hydrogen=true);

	//written binary form. Do not use directly.
	unsigned int to_bin(char* buffer)const;

	//read binary form. Do not use directly.
	void read_bin(char* buffer);

private:
	struct snode {
		int score;
		list<int> children;
		snode() { score = 0; }
		snode(int _score) { score = _score; }
	};							//children tree node.
	vector<snode> atoms;		//atom list
	string symbol;				//string of chain.

	/* functions used in converting */

	//rank atoms by its score, and write into rank
	bool rankby(vector<int>& score, vector<int>& rank);

	//update score by scores of neighbours
	bool renew_score(const int index, const vector<int>& rank, const list<int>& _children_tmp, vector<int>& score);

	//label each atom
	void label(int index, int& crtlabel, vector<list<int> >& bonds, const vector<int>& rank, vector<int>& score);

	//recursive compare function. Do not use directly.
	bool compare(const smiles& sms, int myindex, int hisindex, bool* mark)const;

	/* functions used in compiling */

	//scans the loop of smiles, and write into looproots
	void scanloop(int index, int parent, vector<list<int> >& children_list, vector<pair<int, int> >&, bool* mark)const;

	//recursive compile function. Do not use directly.
	string compilepart(int index, vector<list<int> >& children_list, const vector<string>& loopmark, bool hydrogen)const;

	//interpreter atom symbol
	int read_atom(const string& input, map<int, int>& loops, bool hydrogen=true);

	// turn data into score
	static int to_score(int weight, int hydrogen, int connection);

	// atom score to hydrogen number
	static int score2hydro(int score);
	
	// atom score to weight
	static int score2weight(int score);

};


#endif // !SMILES_H
