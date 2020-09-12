#ifndef CHARMATCHITERATOR_H
#define CHARMATCHITERATOR_H

#include "CharNode.h"

using namespace std;

#include <string>
#include <Rcpp.h>


class CharMatchIterator
{
public:
	CharMatchIterator();
	~CharMatchIterator();

	void Reset();
	bool Next(char ch);
	bool HasAnyNext() const;
	bool IsFullyMatching() const;
	
	void AppendString(const string& str);
	void AppendString(const Rcpp::String& str);
	void AppendStringVector(Rcpp::StringVector& strVec);
	void AppendChars(const char* chStr, size_t chSize);
	
	size_t totalCounter() const;
	void resetTotalCounter();
	void incrementTotalCounter();
	long startMatchIndex() const;
	int matchedCount() const;
private:
	CharNode _root = CharNode('\0');
	CharNode* _currentNode;
	int _matchedChars = 0;
	size_t _counter = 0;
};

#endif // !CHARMATCHITERATOR_H
