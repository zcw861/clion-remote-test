/*
#include "stack_5.1.hpp"

void peek(Stack &st, int index) {
    cout << endl;
    string t;
    if(st.peek(index, t))
        cout << "peek: " << t;
    else
        cout << "peek failed!" << endl;
}

int main() {
    LIFO_Stack lst;
    string str;
    while(cin >> str && !lst.full()) {
        lst.push(str);
    }

    cout << '\n' << "About to call peek() with LIFO_Stack" << endl;
    peek(lst, lst.top()-1);
    cout << lst;

    Peekback_Stack pst;
    while(!lst.empty()) {
        string t;
        if(lst.pop(t))
            pst.push(t);
    }

    cout << '\n' << "About to call peek() with Peekback_Stack" << endl;
    peek(pst, pst.top()-1);
    cout << pst;

    return 0;
}
*/
