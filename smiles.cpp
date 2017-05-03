//a newer version of molecule
#include <stack>
#include <string>
#include <functional>
#include "smiles.h"
#include "tablesort.h"
#include "elementname.h"

#define char2int(x) ((x) - 48)
#define HydrogenName ElementNameL[0]

const int primes[35] = {
	2,3,5,7,11, 13,17,19,23,29, 31,33,37,43,47, 53,59,61,67,71,
	73,79,83,89,97, 101,103,107,109,113, 127,131,137,139,149
};

const char bond_symbol[4] = {
	'~', '-', '=', '@'
};

inline int min(int a, int b) {
	return a < b ? a : b;
}

/* Assistant functions */

inline bool isShort(const char name) {
	return find(ShortNames.begin(), ShortNames.end(), name) != ShortNames.end();
}

inline char getBondSymbol(int order) {
	try {
		return bond_symbol[order];
	}
	catch (const out_of_range&) {
		printf("Warning: bond order %d is not defined", order);
		return '*';
	}
}

// compatibility for old stl
inline string _to_string(int i) {
	return std::to_string(static_cast<long long>(i));
}

atomscore_t smiles::to_score(int weight, int hydrogen, int nhconn, int bond) {
	return hydrogen + (bond + (weight + nhconn*maxatmwht)*maxbond)*maxhnum;
}

int smiles::score2weight(atomscore_t score) {
	return (score / (maxhnum * maxbond)) % maxatmwht;
}
int smiles::score2hydro(atomscore_t score) {
	return score % maxhnum;
}

// interface

void smiles::_push_back(int score) {
	atoms.push_back(snode(score));
}
void smiles::_push_child(int parent, int child, int order) {
	atoms[parent].children.push_back({ child, order });
}

void smiles::push_atom(int weight, int hydrogen) {
	if (hydrogen == -1)hydrogen = DefaultHydrogen.at(weight);
	atoms.push_back(snode(to_score(weight, hydrogen, 0, hydrogen)));
}

void smiles::connect_atom(int atom1, int atom2, int order, bool change_hydrogen) {
	atoms[atom1].children.push_back({ atom2, order });
	atoms[atom2].children.push_back({ atom1, order });
	atoms[atom1].score += maxatmwht * maxhnum * maxbond + maxhnum * order;
	atoms[atom2].score += maxatmwht * maxhnum * maxbond + maxhnum * order;
	if (change_hydrogen){
		atoms[atom1].score -= min(order, score2hydro(atoms[atom1].score));
		atoms[atom2].score -= min(order, score2hydro(atoms[atom2].score));
	}
}

int smiles::element_num(int weight)const {
	int out = 0;
	for (const auto& a : atoms) {
		if (score2weight(a.score) == weight)out++;
	}
	return out;
}

// comparing

bool smiles::operator==(const smiles& sms)const {
	if (atoms.size() != sms.atoms.size())return false;
	bool* mark = new bool[atoms.size()]{ false };
	return compare(sms, 0, 0, mark);
}

bool smiles::compare(const smiles& sms, int myindex, int hisindex, bool* mark)const {
	mark[myindex] = true;
	if (atoms[myindex].score != sms.atoms[hisindex].score)return false;
	if (atoms[myindex].children.size() != sms.atoms[hisindex].children.size())return false;
	for (auto mychild = atoms[myindex].children.begin(), hischild = sms.atoms[hisindex].children.begin();
		mychild != atoms[myindex].children.end(); mychild++, hischild++) {
		if (mark[mychild->first])continue;
		if (mychild->second != hischild->second)return false;
		if (!compare(sms, mychild->first, hischild->first, mark))return false;
	}
	return true;
}

// canonicalize

void smiles::canonicalize()
{
	//1.copy data (score and childtmp is modifyable)
	vector<list<int> > _children_tmp(atoms.size());//neighbour
	vector<atomscore_t> _score(atoms.size()), rank(atoms.size());//score is the unstable score; no matter what happens, score is host. rank is client.
	for (unsigned int i = 0; i < atoms.size(); i++) {
		_children_tmp[i] = atoms[i].neighbour_list();
		_score[i] = atoms[i].score;
	}

	//2.get symmmetry classes (O(N^3))
	while (rankby(_score, rank)) {//if rank is changed
		//O(N*k)
		for (unsigned int i = 0; i < _score.size(); i++) {
			renew_score(i, rank, _children_tmp[i], _score);
		}
	}
	
	//take atom weight into account...
	for (unsigned int i = 0; i < _score.size(); i++) {
		_score[i] = rank[i] * maxatmwht * maxhnum * maxnhcon * maxbond + atoms[i].score;
	}
	rankby(_score, rank);

	//3.canonical labeling
	//find the start
	_score = rank;//here score is served as the symmetry class. (host) rank is served as label. (client)

	int lowest = 0, crtlabel = 1;
	for (const auto& s : _score) {
		if (s > 0)lowest++;
		else break;
	}
	for (auto& r : rank)
		r = -1;
	rank[lowest] = 0;//clear the label vector.

	//tie break algorithm from open babel
	label(lowest, crtlabel, _children_tmp, _score, rank);

	//4.O( Nlog(N) ) rank atoms
	auto rank_cpy = rank;
	sortwith(rank_cpy, atoms, less<atomscore_t>());//put the order of atoms

	//O(N*m) (replace children and rank children, as at this time children is the old index)
	for (auto& atom : atoms) {
		for (auto& child : atom.children)
			child.first = rank[child.first];
		atom.children.sort(greater<pair<int, int> >());
	}
}
bool smiles::rankby(vector<atomscore_t>& score, vector<atomscore_t>& rank) {
	//a more complex algorithm...
	//first we rank the array, then...

	vector<int> sort_result(score.size());
	for (unsigned int i = 0; i < score.size(); i++)
		sort_result[i] = i;
	sortby(sort_result, score, less<atomscore_t>());

	//compare and delete equal rank
	//backup old rank
	vector<atomscore_t> _newrank(rank.size());

	atomscore_t highest_in_old = 0;
	for (const auto& r : rank)
		if (r > highest_in_old)highest_in_old = r;

	_newrank[sort_result[0]] = 0;
	int rank_cnt = 0;
	for (unsigned int i = 1; i < rank.size(); i++) {
		if (score[sort_result[i]] > score[sort_result[i - 1]])rank_cnt++;
		_newrank[sort_result[i]] = rank_cnt;
	}

	if (rank_cnt > highest_in_old) {
		//copy rank
		rank = _newrank; return true;
	}
	else return false;
}
bool smiles::renew_score(const int index, const vector<atomscore_t>& rank, const list<int>& _children_tmp, vector<atomscore_t>& score) {
	atomscore_t oldscore = score[index];
	score[index] = 1;
	for (auto child : _children_tmp) {
		score[index] *= primes[rank[child]];
	}
	return score[index] >= oldscore;
}
void smiles::label(int index, int& crtlabel, vector<list<int> >& children, const vector<atomscore_t>& symclass, vector<atomscore_t>& labels) {
	//remove already labeled
	for (auto child = children[index].begin(); child != children[index].end();) {
		if (labels[*child] >= 0)child = children[index].erase(child);
		else child++;
	}
	sortby(children[index], symclass, greater<atomscore_t>());

	//then, label children
	for (const auto& child : children[index]) {
		labels[child] = crtlabel;
		crtlabel++;
	}

	//finally, label others
	for (const auto& child : children[index])
		label(child, crtlabel, children, symclass, labels);
}


// io functions
string smiles::to_smiles(bool hydrogen) {
	if (this->symbol.empty())
		compile(hydrogen);
	return symbol;
}

void smiles::compile(bool hydrogen, bool single_bond)
{
	//first dfs: scan loop&delete parent
	vector<pair<int, int> > loops;
	vector<list<pair<int, int> > > children_list(atoms.size(), list<pair<int, int> >());
	bool* mark = new bool[atoms.size()]{ 0 };

	scanloop(0, -1, children_list, loops, mark);
	delete[] mark;

	vector<string> loopmark(atoms.size());
	int crtloopid = 1;
	for (unsigned int i = 0; i < loops.size(); i++) {
		if (!loopmark[loops[i].second].empty()) {
			loopmark[loops[i].second].push_back(',');
		}
		loopmark[loops[i].second] += _to_string(i + 1);
	}

	symbol = compilepart(0, children_list, loopmark, hydrogen, single_bond);
}
void smiles::scanloop(int index, int parent, vector<list<pair<int, int> > >& children_list, vector<pair<int, int> >& loops, bool* mark)const
{
	mark[index] = true;
	for (const auto& child : atoms[index].children) {
		if (child.first == parent) {
			continue;
		}
		else {
			if (mark[child.first]) {
				if (find(loops.begin(), loops.end(), pair<int, int>{child.first, index}) == loops.end()) {
					loops.push_back(pair<int, int>{index, child.first});
					children_list[index].push_back({ -(int)loops.size(), child.second });	//use minus number to represent loop
				}
			}
			else {
				children_list[index].push_back(child);
				scanloop(child.first, index, children_list, loops, mark);
			}
		}
	}
}
string smiles::compilepart(int index, vector<list<pair<int, int> > >& children_list, const vector<string>& loopmark, bool hydrogen, bool single_bond)const {
	//roots
	if (index < 0) {
		return _to_string(-index);
	}

	string result = "";
	result.append(WeightNameL.at(score2weight(atoms[index].score)));
	result.append(loopmark[index]);

	if (hydrogen) {
		string hydrogen_str;
		if (score2hydro(atoms[index].score) > 1)hydrogen_str = HydrogenName + _to_string(score2hydro(atoms[index].score));
		else if (score2hydro(atoms[index].score) == 1)hydrogen_str = HydrogenName;

		if (!hydrogen_str.empty()) {
			result = "[" + result + hydrogen_str + "]";
		}
	}

	for (auto child = children_list[index].rbegin(); child != children_list[index].rend(); child++) {
		if (child == --children_list[index].rend()) {
			if (single_bond || child->second != 1) {
				result.push_back(getBondSymbol(child->second));
			}
			result += compilepart(child->first, children_list, loopmark, hydrogen, single_bond);
		}
		else {
			result += "(";
			if (single_bond || child->second != 1) {
				result.push_back(getBondSymbol(child->second));
			}
			result += compilepart(child->first, children_list, loopmark, hydrogen, single_bond) + ")";
		}
	}
	return result;
}

void smiles::read_smiles(const string& input, bool read_hydrogen) {
	stack<string::const_iterator> braket_stack;
	stack<int> atom_stack;
	map<int, int> loop_roots;

	for (auto s = input.begin(); s != input.end(); s++) {
		if (*s == '[')braket_stack.push(s);
		else if (*s == ']') {
			auto last_bra_pos = braket_stack.top();
			braket_stack.pop();
			if (braket_stack.empty()) {
				atoms.push_back(snode(read_atom(string(last_bra_pos, s + 1), loop_roots, read_hydrogen)));
				if (!atom_stack.empty()) {
					connect_atom(atom_stack.top(), atoms.size() - 1, 1, false);
					atom_stack.pop();
				}
				atom_stack.push(atoms.size() - 1);
			}
		}
		else if (*s == ')') {
			atom_stack.pop();
		}
		else if (braket_stack.empty() && isShort(*s)) { // trivial
			atoms.push_back(snode(read_atom(string(s, s + 1), loop_roots, read_hydrogen)));
			if (!atom_stack.empty()) {
				connect_atom(atom_stack.top(), atoms.size() - 1, 1, false);
				atom_stack.pop();
			}			
			atom_stack.push(atoms.size() - 1);
		}
		else if (braket_stack.empty() && *s >= '0' && *s <= '9') {
			auto root = loop_roots.find(char2int(*s));
			if (root == loop_roots.end())loop_roots[char2int(*s)] = atoms.size() - 1;
			else {
				connect_atom(root->second, atoms.size() - 1, 1, false);
			}
		}
	}
}

atomscore_t smiles::read_atom(const string& input, map<int, int>& loops, bool read_hydrogen) {
	string atom_str;
	string::const_iterator atom_end;

	if (input[0] == '[') {
		if (input[1] == '[') {
			atom_end = ++find(input.begin(), input.end(), ']');
			atom_str = string(input.begin() + 1, atom_end);
		}
		else if (isShort(input[1])) {
			atom_end = input.begin() + 2;
			atom_str = string(1, input[1]);
		}
		else {
			atom_end = ++find(input.begin(), input.end(), ']');
			atom_str = string(input.begin(), atom_end);
		}
	}
	else {
		atom_end = input.begin() + 1;
		atom_str = string(1, input[0]);
	}

	if (atom_end != input.end() && *atom_end >= '0' && *atom_end <= '9') {
		loops[char2int(*atom_end)] = atoms.size();
		atom_end++;
	}

	int weight = NameLWeight.at(atom_str);
	int hydrogen = 0;
	if (!read_hydrogen) {
		hydrogen = DefaultHydrogen.at(weight);
	}
	else if (atom_end != input.end() && *atom_end == 'H'){
		if (*(atom_end + 1) == ']')hydrogen = 1;
		else hydrogen = char2int(*(atom_end + 1));
	}
	return to_score(weight, hydrogen, 0, hydrogen);
}

list<int> smiles::snode::neighbour_list() const
{
	list<int> neighbours;
	for (const auto& child : children) {
		neighbours.push_back(child.first);
	}
	return neighbours;
}
