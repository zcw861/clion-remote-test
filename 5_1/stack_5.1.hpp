//
// Created by Administrator on 2025/8/3.
//

#ifndef STACK_5.1_HPP
#define STACK_5.1_HPP

#include <string>
#include <iostream>
#include <vector>

using namespace std;

typedef string elemtype;

class Stack {
public:
    virtual ~Stack() {}

    virtual bool empty() const = 0;
    virtual bool full() const = 0;

    virtual int size() const = 0;

    virtual bool pop(elemtype&) = 0;
    virtual bool push(elemtype&) = 0;
    virtual bool peek(int index, elemtype&) const = 0;

    virtual void print(ostream& = cout) const = 0;
};

ostream& operator<< (ostream &os, const Stack &rhs) {rhs.print(); return os;}

class LIFO_Stack : public Stack {
public:
    LIFO_Stack(int capacity = 0) : _top(0) {
        if(capacity)
            _stack.reserve(capacity);
    }

    int size() const {return _stack.size();}
    bool empty() const {return _stack.empty();}
    bool full() const {return size() >= _stack.max_size();}
    int top() const {return _top;}

    void print(ostream &os = cout) const;

    bool pop(elemtype &elem);
    bool push(elemtype &elem);
    bool peek(int, elemtype&) const {return false;}

private:
    int _top;
    vector<elemtype> _stack;
};

inline bool LIFO_Stack::pop(elemtype &elem) {
    if(empty()) return false;
    elem = _stack[--_top];
    _stack.pop_back();
    return true;
}

inline bool LIFO_Stack::push(elemtype &elem) {
    if(full()) return false;
    _stack.push_back(elem);
    ++_top;
    return true;
}

inline void LIFO_Stack::print(ostream &os) const {
    vector<elemtype>::const_reverse_iterator
        rit = _stack.rbegin(),
        rend = _stack.rend();

    os << "\n\t";
    while (rit != rend) {
        os << *rit++ << "\n\t";
    }
    os << endl;
}

class Peekback_Stack : public Stack {
public:
    Peekback_Stack(int capacity = 0) : _top(0) {
        if(capacity)
            _stack.reserve(capacity);
    }

    int top() const {return _top;}
    int size() const {return _stack.size();}
    bool empty() const {return _stack.empty();}
    bool full() const {return size() >= _stack.max_size();}

    void print(ostream &os = cout) const;

    bool pop(elemtype &elem);
    bool push(elemtype &elem);
    bool peek(int index, elemtype &elem) const;

private:
    int _top;
    vector<elemtype> _stack;
};

inline bool Peekback_Stack::pop(elemtype &elem) {
    if(empty()) return false;
    elem = _stack[--_top];
    _stack.pop_back();
    return true;
}

inline bool Peekback_Stack::push(elemtype &elem) {
    if(full()) return false;
    _stack.push_back(elem);
    ++_top;
    return true;
}

inline bool Peekback_Stack::peek(int index, elemtype &elem) const {
    if(empty() || index < 0 || index >= size()) return false;
    elem = _stack[index];
    return true;
}

inline void Peekback_Stack::print(ostream &os) const {
    vector<elemtype>::const_reverse_iterator
        rit = _stack.rbegin(),
        rend = _stack.rend();

    os << "\n\t";
    while (rit != rend) {
        os << *rit++ << "\n\t";
    }
    os << endl;
}

#endif //STACK_5.1_HPP
