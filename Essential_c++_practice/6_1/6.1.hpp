//
// Created by Administrator on 2025/9/13.
//

#ifndef INC_6_1_HPP
#define INC_6_1_HPP

template<typename elemType>
class example {
public:
    example(const example &min, const example &max);
    example(const example* array, int size);

    elemType& operator[](int index);
    bool operator==(const example&) const;

    bool insert(const elemType*, int);
    bool insert(const elemType&);

    elemType min() const {return _min;}
    elemType max() const {return _max;}

    void min(const elemType&);
    void max(const elemType&);

    int count(const elemType &value) const;

private:
    int _size;
    elemType *_parray;
    elemType _min;
    elemType _max;
};


#endif //INC_6_1_HPP
