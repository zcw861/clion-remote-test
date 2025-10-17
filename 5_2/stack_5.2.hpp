//
// Created by Administrator on 2025/8/4.
//

#ifndef STACK_5.2_HPP
#define STACK_5.2_HPP

#include <string>
#include <iostream>
#include <vector>

using namespace std;

typedef string elemtype;

class Stack {
public:
    Stack(int capacity = 0) : _top(0) {
        if(capacity)
            _stack.reserve(capacity);
    }

    void print(ostream &os = cout) const;

    virtual ~Stack() {}

    int size() const {return _stack.size();}
    bool empty() const {return _stack.empty();}
    bool full() const {return size() >= _stack.max_size();}
    int top() const {return _top;}

    bool pop(elemtype &elem);
    bool push(elemtype &elem);
    virtual bool peek(int, elemtype&) const {return false;}

protected:
    int _top;
    vector<elemtype> _stack;
};

ostream& operator<< (ostream &os, const Stack &rhs) {rhs.print(); return os;}

inline bool Stack::pop(elemtype &elem) {
    if(empty()) return false;
    elem = _stack[--_top];
    _stack.pop_back();
    return true;
}

inline bool Stack::push(elemtype &elem) {
    if(full()) return false;
    _stack.push_back(elem);
    ++_top;
    return true;
}

inline void Stack::print(ostream &os) const {
    vector<elemtype>::const_reverse_iterator
        rit = _stack.rbegin(),
        rend = _stack.rend();

    os << "\n\t";
    while (rit != rend) {
        os << *rit++ << "\n\t";
    }
    os << endl;
}


class LIFO_Stack : public Stack {
public:
    LIFO_Stack(int capacity = 0) : Stack(capacity) {}

    virtual bool peek(int, elemtype&) const {return false;}
};



class Peekback_Stack : public Stack {
public:
    Peekback_Stack(int capacity = 0) : Stack(capacity) {}

    virtual bool peek(int index, elemtype &elem) const;
};

inline bool Peekback_Stack::peek(int index, elemtype &elem) const {
    if(empty() || index < 0 || index >= size()) return false;
    elem = _stack[index];
    return true;
}

#endif //STACK_5.2_HPP
