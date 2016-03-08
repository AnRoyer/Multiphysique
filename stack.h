#ifndef STACK_H_INCLUDED
#define STACK_H_INCLUDED

#include <string>

struct Cell
{
    struct Cell *after;
    std::string name;
};
typedef struct Cell Cell;

class Stack
{
private:

    Cell *up;

public:

    Stack();
    void push(std::string);
    std::string peek();
    std::string pop();
    bool isEmpty() const;
    void clear();
    ~Stack();

};

#endif // STACK_H_INCLUDED
