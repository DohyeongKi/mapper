#ifndef PAIR_H_
#define PAIR_H_

class Pair {
public:
	int data;
	char type;

	Pair(int it, char tp);
	Pair(int it);
	Pair();
};

Pair::Pair(int it, char tp){
	data = it;
	type = tp;
}

Pair::Pair(int it){
	data = it;
	type = 'A';
}

Pair::Pair(){
	data = -1;
	type = 'A';
}

#endif /* PAIR_H_ */

/*
 Type - A: not determined, B: only that level
        C: that + upper level, D: that + lower level
        E: that + upper + lower level
 */
