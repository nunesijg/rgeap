#include "CharNode.h"

CharNode::CharNode(const char ch)
{
	_char = ch;
	_endOfString = false;
}

CharNode::CharNode(const char ch, const bool endOfString)
{
	_char = ch;
	_endOfString = endOfString;
}


CharNode::~CharNode()
{
	for (auto it = _nextNodes.begin(); it != _nextNodes.end(); ++it)
	{
		delete it->second;
	}
	_nextNodes.clear();
}

CharNode* CharNode::appendChar(const char ch, const bool end)
{
	CharNode* ret;
	const auto it = _nextNodes.find(ch);
	if (it == _nextNodes.end())
	{
		ret = new CharNode(ch, end);
		_nextNodes[ch] = ret;
	}
	else
	{
		ret = it->second;
		if (end)
			ret->_endOfString = true;
	}
	return ret;
}

CharNode* CharNode::getNext(const char ch)
{
	const auto it = _nextNodes.find(ch);
	if (it == _nextNodes.end()) return nullptr;
	return it->second;
}

char CharNode::getChar() const
{
	return _char;
}

bool CharNode::anyNext() const
{
	return !_nextNodes.empty();
}

bool CharNode::isEndOfString() const
{
	return _endOfString;
}
