#ifndef CHARNODE_H
#define CHARNODE_H

using namespace std;

#include <map>

class CharNode
{
public:
	explicit CharNode(char ch);
	explicit CharNode(char ch, bool endOfString);
	~CharNode();

	CharNode* appendChar(char ch, bool end);
	CharNode* getNext(char ch);
	char getChar() const;
	bool anyNext() const;
	bool isEndOfString() const;
private:
	map<char, CharNode*> _nextNodes;
	char _char;
	bool _endOfString = false;
};

#endif // !CHARNODE_H
