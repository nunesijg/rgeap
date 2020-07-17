#include "CharMatchIterator.h"
#include "funhelpers.h"

CharMatchIterator::CharMatchIterator()
{
	_currentNode = &_root;
	_matchedChars = 0;
}


CharMatchIterator::~CharMatchIterator()
{
}

void CharMatchIterator::Reset()
{
	_currentNode = &_root;
	_matchedChars = 0;
}

bool CharMatchIterator::Next(const char ch)
{
	const auto next = _currentNode->getNext(ch);
	_counter++;
	if (next == nullptr) return false;
	_currentNode = next;
	_matchedChars++;
	return true;
}

bool CharMatchIterator::HasAnyNext() const
{
	return _currentNode->anyNext();
}

bool CharMatchIterator::IsFullyMatching() const
{
	return _currentNode->isEndOfString();
}

void CharMatchIterator::AppendString(const string& str)
{
  const char* chStr = str.c_str();
  const size_t chSize = str.size();
  AppendChars(chStr, chSize);
}

void CharMatchIterator::AppendString(const Rcpp::String& str)
{
  if (str == NA_STRING) return;
  const size_t chSize = str_size(str);
  if (chSize == 0) return;
  const char* chStr = str.get_cstring();
  AppendChars(chStr, chSize);
}

void CharMatchIterator::AppendStringVector(Rcpp::StringVector& strVec)
{
  Rcpp::StringVector::iterator it = strVec.begin();
  while (it != strVec.end())
  {
    Rcpp::String str = *it;
    AppendString(str);
    ++it;
  }
}

void CharMatchIterator::AppendChars(const char* chStr, const size_t chSize)
{
  CharNode* node = &_root;
  for (size_t i = 0; i < chSize; i++)
  {
    const bool strEnd = i == chSize - 1;
    node = node->appendChar(chStr[i], strEnd);
  }
}

size_t CharMatchIterator::totalCounter() const
{
	return _counter;
}

void CharMatchIterator::incrementTotalCounter()
{
  _counter++;
}

void CharMatchIterator::resetTotalCounter()
{
	_counter = 0;
}

long CharMatchIterator::startMatchIndex() const
{
	return static_cast<long>(_counter) - _matchedChars;
}

int CharMatchIterator::matchedCount() const
{
	return _matchedChars;
}

