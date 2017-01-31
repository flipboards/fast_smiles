#pragma once
#include<vector>
#include<list>

using namespace std;

//sort a index vector by score using bubble sort. Pred means the order of head to tail.
template<typename type, typename order>
void sortby(vector<type>& source, const vector<int>& score, const order pred) {
	for (unsigned int i = 0; i < source.size() - 1; i++) {
		for (auto source_j = source.begin(); source_j < source.end() - 1 - i; source_j++) {
			if (pred(score[*(source_j+1)], score[*source_j])) {
				swap(*source_j, *(source_j + 1));
			}
		}
	}
}

//sort a index list by score using bubble sort. Pred means the order of head to tail.
template<typename type, typename order>
void sortby(list<type>& source, const vector<int>& score, const order pred) {
	if (source.empty())return;
	auto sort_end = source.end();
	while(sort_end != source.begin()) {
		auto source_j = source.begin();
		auto next_j = source_j;
		next_j++;

		while(next_j != sort_end) {
			if (pred(score[*next_j], score[*source_j])) {
				swap(*source_j, *next_j);
			}
			source_j = next_j;
			next_j++;
		}
		sort_end--;
	}
}

//rank by a exist rank
template<typename iter>
void sortwithrank(iter begin, iter end, vector<int>& rank) {
	if (end == begin || end == begin + 1)return;
	unsigned int i = 0;
	for (iter s = begin; s != end; s++) {
		if (rank[i] != i) {
			swap(*s, *(begin + rank[i]));
			int t = rank[rank[i]];
			rank[rank[i]] = rank[i];
			rank[i] = t;
		}
		i++;
	}
}

//sort a list with another list accompanied. Pred means the order of head to tail.
template<typename type1, typename type2, typename order>
void sortwith(vector<type1>& score, vector<type2>& assist, const order pred) {
	for (size_t i = 0; i < score.size() - 1; i++) {
		for (size_t j = 0; j < score.size() - 1 - i; j++) {
			if (pred(score[j + 1], score[j])) {
				swap(score[j], score[j + 1]);
				swap(assist[j], assist[j + 1]);
			}
		}
	}

}
