//a newer version of molecule
#include"smiles.h"
#include<stack>
#include<string>
#include<functional>
#include"tablesort.h"
#include"elementname.h"

// #define lettermode	1
// #define hmode		2

#define char2int(x) ((x) - 48)
#define HydrogenName ElementNameL[0]

const int primes[35] = {
	2,3,5,7,11, 13,17,19,23,29, 31,33,37,43,47, 53,59,61,67,71,
	73,79,83,89,97, 101,103,107,109,113, 127,131,137,139,149
};

/* Assistant functions */

inline bool isShort(const char name) {
	return find(ShortNames.begin(), ShortNames.end(), name) != ShortNames.end();
}

int to_score(int weight, int hydrogen, int connection) {
	return hydrogen + (weight + connection*maxatmwht)*maxhnum;
}

inline int score2weight(int score) {
	return (score / maxhnum) % maxatmwht;
}
inline int score2hydro(int score) {
	return score % maxhnum;
}

// interface

void smiles::push_back(int score) {
	atoms.push_back(snode(score));
}

void smiles::push_atom(int weight, int hydrogen) {
	if (hydrogen == -1)hydrogen = DefaultHydrogen.at(weight);
	atoms.push_back(snode(to_score(weight, hydrogen, 0)));
}

void smiles::connect_atom(int atom1, int atom2, bool change_hydrogen) {
	atoms[atom1].children.push_back(atom2);
	atoms[atom2].children.push_back(atom1);
	atoms[atom1].score += maxatmwht * maxhnum;
	atoms[atom2].score += maxatmwht * maxhnum;
	if (change_hydrogen){
		if (score2hydro(atoms[atom1].score) > 0)atoms[atom1].score--;
		if (score2hydro(atoms[atom2].score) > 0)atoms[atom2].score--;
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
		if (mark[*mychild])continue;
		if (!compare(sms, *mychild, *hischild, mark))return false;
	}
	return true;
}

// canonicalize

void smiles::canonicalize()
{
	//1.copy data (score and childtmp is modifyable)
	vector<list<int> > _children_tmp(atoms.size());//neighbour
	vector<int> _score(atoms.size()), rank(atoms.size());//score is the unstable score; no matter what happens, score is host. rank is client.
	for (unsigned int i = 0; i < atoms.size(); i++) {
		_children_tmp[i] = atoms[i].children;
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
		_score[i] = rank[i] * maxatmwht * maxhnum * maxnhcon + atoms[i].score;
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
	sortwith(rank_cpy, atoms, less<int>());//put the order of atoms

	//O(N*m) (replace children and rank children, as at this time children is the old index)
	for (auto& atom : atoms) {
		for (auto& child : atom.children)
			child = rank[child];
		atom.children.sort(greater<int>());
	}
}
bool smiles::rankby(vector<int>& score, vector<int>& rank) {
	//a more complex algorithm...
	//first we rank the array, then...

	vector<int> sort_result(score.size());
	for (unsigned int i = 0; i < score.size(); i++)
		sort_result[i] = i;
	sortby(sort_result, score, less<int>());

	//compare and delete equal rank
	//backup old rank
	vector<int> _newrank(rank.size());

	int highest_in_old = 0;
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
bool smiles::renew_score(const int index, const vector<int>& rank, const list<int>& _children_tmp, vector<int>& score) {
	int oldscore = score[index];
	score[index] = 1;
	for (auto child : _children_tmp) {
		score[index] *= primes[rank[child]];
	}
	return score[index] >= oldscore;
}
void smiles::label(int index, int& crtlabel, vector<list<int> >& children, const vector<int>& symclass, vector<int>& labels) {
	//remove already labeled
	for (auto child = children[index].begin(); child != children[index].end();) {
		if (labels[*child] >= 0)child = children[index].erase(child);
		else child++;
	}
	sortby(children[index], symclass, greater<int>());

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

void smiles::compile(bool hydrogen)
{
	//first dfs: scan loop&delete parent
	vector<pair<int, int> > loops;
	vector<list<int> > children_list(atoms.size(), list<int>());
	bool* mark = new bool[atoms.size()]{ 0 };

	scanloop(0, -1, children_list, loops, mark);
	delete[] mark;

	vector<string> loopmark(atoms.size());
	int crtloopid = 1;
	for (unsigned int i = 0; i < loops.size(); i++)
		loopmark[loops[i].second] += std::to_string(i + 1);

	symbol = compilepart(0, children_list, loopmark, hydrogen);
}
void smiles::scanloop(int index, int parent, vector<list<int> >& children_list, vector<pair<int, int> >& loops, bool* mark)const
{
	mark[index] = true;
	for (const auto& child : atoms[index].children) {
		if (child == parent) {
			continue;
		}
		else {
			if (mark[child]) {
				if (find(loops.begin(), loops.end(), pair<int, int>{child, index}) == loops.end()) {
					loops.push_back(pair<int, int>{index, child});
					children_list[index].push_back(-(int)loops.size());	//use minus number to represent loop
				}
			}
			else {
				children_list[index].push_back(child);
				scanloop(child, index, children_list, loops, mark);
			}
		}
	}
}
string smiles::compilepart(int index, vector<list<int> >& children_list, const vector<string>& loopmark, bool hydrogen)const {
	//roots
	if (index < 0) {
		return std::to_string(-index);
	}

	string result = "";
	result.append(WeightNameL.at(score2weight(atoms[index].score)));
	result.append(loopmark[index]);

	if (hydrogen) {
		string hydrogen_str;
		if (score2hydro(atoms[index].score) > 1)hydrogen_str = HydrogenName + std::to_string(score2hydro(atoms[index].score));
		else if (score2hydro(atoms[index].score) == 1)hydrogen_str = HydrogenName;

		if (!hydrogen_str.empty()) {
			result = "[" + result + hydrogen_str + "]";
		}
	}

	for (auto child = children_list[index].rbegin(); child != children_list[index].rend(); child++) {
		if (child == --children_list[index].rend())
			result += compilepart(*child, children_list, loopmark, hydrogen);
		else
			result += "(" + compilepart(*child, children_list, loopmark, hydrogen) + ")";
	}
	return result;
}
/*
void smiles::read_smiles(const string& symbol) {
	//first read the atom, then scan the connectivity

	struct bond { int id_1, id_2; };
	vector<bond> loops;

	enum{nomode, lettermode, hmode};
	auto mode = lettermode;

	stack<size_t> atomstack;
	string name_buffer;
	int prev_cnt = -1;

	auto s = symbol.begin();
	while (s != symbol.end()) {
		if (*s >= 'A' && *s != 'H' && mode == lettermode)name_buffer.push_back(*s);
		else if (*s < 'A' && *s >= '0') {
			//end letter
			switch (mode)
			{
			case lettermode:
				mode = nomode;
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
			case nomode:
				if (loops.size() < char2int(*s)) {
					loops.resize(char2int(*s)); loops[char2int(*s)].id_1 = atoms.size() - 1;
				}
				else loops[char2int(*s)].id_2 = atoms.size() - 1;
				break;
			case hmode:
				atoms.back().score += char2int(*s);
				mode = nomode;
				break;
			default:
				break;
			}
		}
		else if (*s == 'h') {
			if (mode == lettermode) {
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
			}
			mode = hmode;
		}
		else if (*s == '-' || *s == '(' || *s == ')') {
			if (mode == lettermode) {
				atomstack.push(atoms.size());
				atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
				name_buffer.clear();
				if (prev_cnt >= 0) {
					atoms.back().children.push_back(prev_cnt);
					atoms[prev_cnt].children.push_back(atoms.size() - 1);
				}
			}
			else if (mode == hmode) {
				atoms.back().score += 1;
				mode = nomode;
			}
			prev_cnt = atomstack.top(); mode = lettermode;
			if (*s != '(')atomstack.pop();
		}
		s++;
	}
	if (mode == lettermode) {
		atoms.push_back(snode(maxhnum * nameL2weight(name_buffer)));
		name_buffer.clear();
		if (prev_cnt >= 0) {
			atoms.back().children.push_back(prev_cnt);
			atoms[prev_cnt].children.push_back(atoms.size() - 1);
		}
	}
	else if (mode == hmode) {
		atoms.back().score += 1;
		mode = nomode;
	}

	for (auto loop : loops) {
		atoms[loop.id_1].children.push_back(loop.id_2);
		atoms[loop.id_2].children.push_back(loop.id_1);
	}
	for (auto atom = atoms.begin(); atom != atoms.end(); atom++) {
		atom->score += maxhnum*maxatmwht*atom->children.size();
	}
}*/

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
					connect_atom(atom_stack.top(), atoms.size() - 1, false);
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
				connect_atom(atom_stack.top(), atoms.size() - 1, false);
				atom_stack.pop();
			}			
			atom_stack.push(atoms.size() - 1);
		}
		else if (braket_stack.empty() && *s >= '0' && *s <= '9') {
			auto root = loop_roots.find(char2int(*s));
			if (root == loop_roots.end())loop_roots[char2int(*s)] = atoms.size() - 1;
			else {
				connect_atom(root->second, atoms.size() - 1, false);
			}
		}
	}
}

int smiles::read_atom(const string& input, map<int, int>& loops, bool read_hydrogen) {
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
	return to_score(weight, hydrogen, 0);
}

// binary io functions

unsigned int smiles::to_bin(char* buffer)const {
	int* p = (int*)buffer;
	*p = atoms.size(); p++;
	for (const auto& atom : atoms) {
		*p = atom.score; p++;
		*p = atom.children.size(); p++;
		for (const auto& child : atom.children) {
			*p = child; p++;
		}
	}
	return (p - (int*)buffer) * sizeof(int);
}
void smiles::read_bin(char* buffer) {
	int* p = (int*)buffer;
	atoms.resize(*p); p++;
	for (auto& atom : atoms) {
		atom.score = *p; p++;
		atom.children.resize(*p); p++;
		for (auto& child : atom.children) {
			child = *p; p++;
		}
	}
}

